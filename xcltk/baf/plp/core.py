# core.py

import math
import os
import pickle
import pysam
import sys

from .mcount import MCount
from .sam import check_read, sam_fetch
from .zfile import zopen, ZF_F_GZIP

def sp_region(reg, conf):
    reg_ref_umi = {smp:set() for smp in conf.barcodes}
    reg_alt_umi = {smp:set() for smp in conf.barcodes}
    reg_oth_umi = {smp:set() for smp in conf.barcodes}
    mcnt = MCount(conf.barcodes, conf)

    for snp in reg.snp_list:
        itr = sam_fetch(conf.sam, snp.chrom, snp.pos, snp.pos)
        if not itr:    
            continue
        if mcnt.add_snp(snp) < 0:   # mcnt reset() inside.
            return((-3, None, None, None, None))
        for read in itr:
            if check_read(read, conf) < 0:
                continue
            ret = mcnt.push_read(read)
            if ret < 0:
                if ret == -1:
                    return((-5, None, None, None, None))
                continue
        if mcnt.stat() < 0:
            return((-7, None, None, None, None))
        snp_cnt = sum(mcnt.tcount)
        if snp_cnt < conf.min_count:
            continue
        snp_ref_cnt = mcnt.tcount[mcnt.base_idx[snp.ref]]
        snp_alt_cnt = mcnt.tcount[mcnt.base_idx[snp.alt]]
        snp_minor_cnt = min(snp_ref_cnt, snp_alt_cnt)
        if snp_minor_cnt < snp_cnt * conf.min_maf:
            continue
        for smp, scnt in mcnt.cell_cnt.items():
            for umi, ucnt in scnt.umi_cnt.items():
                if not ucnt.allele:
                    continue
                ale_idx = snp.get_region_allele_index(ucnt.allele)
                if ale_idx == 0:        # ref allele of the region.
                    reg_ref_umi[smp].add(umi)
                elif ale_idx == 1:      # alt allele of the region.
                    reg_alt_umi[smp].add(umi)
                else:
                    reg_oth_umi[smp].add(umi)

    reg_ref_cnt = {smp:0 for smp in conf.barcodes}
    reg_alt_cnt = {smp:0 for smp in conf.barcodes}
    reg_oth_cnt = {smp:0 for smp in conf.barcodes}
    reg_dp_cnt =  {smp:0 for smp in conf.barcodes}
    for smp in conf.barcodes:
        reg_ref_cnt[smp] = len(reg_ref_umi[smp])
        reg_alt_cnt[smp] = len(reg_alt_umi[smp])
        dp_umi = reg_ref_umi[smp].union(reg_alt_umi[smp]) # CHECK ME! theoretically no shared UMIs
        reg_oth_umi[smp] = reg_oth_umi[smp].difference(dp_umi)
        reg_oth_cnt[smp] = len(reg_oth_umi[smp])
        reg_dp_cnt[smp]  = len(dp_umi)
    
    return((0, reg_ref_cnt, reg_alt_cnt, reg_oth_cnt, reg_dp_cnt, reg_ref_umi, reg_alt_umi, reg_oth_umi))

# TODO: use clever IPC (Inter-process communication) instead of naive `raise Error`.
# NOTE: 
# 1. bgzf errors when using pysam.AlignmentFile.fetch in parallel (with multiprocessing)
#    https://github.com/pysam-developers/pysam/issues/397
def sp_count(thdata):
    func = "sp_count"
    conf = thdata.conf
    thdata.ret = -1

    conf.sam = pysam.AlignmentFile(conf.sam_fn, "r")    # auto detect file format

    reg_list = None
    if thdata.is_reg_pickle:
        with open(thdata.reg_obj, "rb") as fp:
            reg_list = pickle.load(fp)
        os.remove(thdata.reg_obj)
    else:
        reg_list = thdata.reg_obj

    fp_reg = zopen(thdata.out_region_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_ad = zopen(thdata.out_ad_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_dp = zopen(thdata.out_dp_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_oth = zopen(thdata.out_oth_fn, "wt", ZF_F_GZIP, is_bytes = False)

    m_reg = float(len(reg_list))
    n_reg = 0
    l_reg = 0
    k_reg = 1
    for reg_idx, reg in enumerate(reg_list):
        if conf.debug > 0:
            sys.stderr.write("[D::%s][Thread-%d] processing region '%s' ...\n" %
                              (func, thdata.idx, reg.name))

        if reg.snp_list:
            ret, reg_ref_cnt, reg_alt_cnt, reg_oth_cnt, reg_dp_cnt, \
                 reg_ref_umi, reg_alt_umi, reg_oth_umi = sp_region(reg, conf)
            if ret < 0:
                raise ValueError("[%s] errcode %d" % (func, -9))

            str_reg, str_ad, str_dp, str_oth = "", "", "", ""
            for i, smp in enumerate(conf.barcodes):
                nu_ad, nu_dp, nu_oth = -1, -1, -1
                if reg_ref_cnt[smp] + reg_alt_cnt[smp] != reg_dp_cnt[smp]:
                    if conf.debug > 0:
                        msg = "[D::%s][Thread-%d] region '%s', sample '%s':\n" % (
                                func, thdata.idx, reg.name, smp)
                        msg += "\tduplicate UMIs: REF, ALT, DP_uniq (%d, %d, %d)!\n" % (
                                reg_ref_cnt[smp], reg_alt_cnt[smp], reg_dp_cnt[smp])
                        sys.stderr.write(msg)
                    if conf.no_dup_hap:
                        nu_share = reg_ref_cnt[smp] + reg_alt_cnt[smp] - reg_dp_cnt[smp]
                        nu_ad = reg_alt_cnt[smp] - nu_share
                        nu_dp = reg_dp_cnt[smp] - nu_share
                        umi_share = reg_ref_umi[smp].intersection(reg_alt_umi[smp])
                        reg_ref_umi[smp] = reg_ref_umi[smp].difference(umi_share)
                        reg_alt_umi[smp] = reg_alt_umi[smp].difference(umi_share)
                    else:
                        nu_ad = reg_alt_cnt[smp]
                        nu_dp = reg_ref_cnt[smp] + reg_alt_cnt[smp]
                else:
                    nu_ad, nu_dp = reg_alt_cnt[smp], reg_dp_cnt[smp]
                nu_oth = reg_oth_cnt[smp]

                if nu_dp + nu_oth <= 0:
                    continue
                if nu_ad > 0:
                    str_ad += "%d\t%d\t%d\n" % (k_reg, i + 1, nu_ad)
                    thdata.nr_ad += 1
                if nu_dp > 0:
                    str_dp += "%d\t%d\t%d\n" % (k_reg, i + 1, nu_dp)
                    thdata.nr_dp += 1
                if nu_oth > 0:
                    str_oth += "%d\t%d\t%d\n" % (k_reg, i + 1, nu_oth)
                    thdata.nr_oth += 1

            if str_dp or str_oth:
                fp_ad.write(str_ad)
                fp_dp.write(str_dp)
                fp_oth.write(str_oth)
                fp_reg.write("%s\t%d\t%d\t%s\n" % (reg.chrom, reg.start, reg.end - 1, reg.name))
                k_reg += 1
            elif conf.output_all_reg:
                fp_reg.write("%s\t%d\t%d\t%s\n" % (reg.chrom, reg.start, reg.end - 1, reg.name))
                k_reg += 1

            if str_dp:
                out_umi_fn = os.path.join(conf.umi_dir, "%s.tsv" % reg.name)
                with open(out_umi_fn, "w") as fp_umi:
                    for smp in conf.barcodes:
                        if smp in reg_ref_umi:
                            umi_set = reg_ref_umi[smp]
                            if umi_set:
                                for umi in umi_set:
                                    fp_umi.write("%s\t%s\t%d\n" % (smp, umi, 0))
                        if smp in reg_alt_umi:
                            umi_set = reg_alt_umi[smp]
                            if umi_set:
                                for umi in umi_set:
                                    fp_umi.write("%s\t%s\t%d\n" % (smp, umi, 1))

        elif conf.output_all_reg:
            fp_reg.write("%s\t%d\t%d\t%s\n" % (reg.chrom, reg.start, reg.end - 1, reg.name))
            k_reg += 1

        n_reg += 1
        frac_reg = n_reg / m_reg
        if frac_reg - l_reg >= 0.02 or n_reg == m_reg:
            sys.stdout.write("[I::%s][Thread-%d] %d%% genes processed\n" % 
                (func, thdata.idx, math.floor(frac_reg * 100)))
            l_reg = frac_reg

    thdata.nr_reg = k_reg - 1

    fp_reg.close()
    fp_ad.close()
    fp_dp.close()
    fp_oth.close()
    conf.sam.close()

    thdata.conf = None    # sam object cannot be pickled.
    thdata.ret = 0

    if thdata.out_fn:
        with open(thdata.out_fn, "wb") as fp_td:
            pickle.dump(thdata, fp_td)
            
    return((0, thdata))
