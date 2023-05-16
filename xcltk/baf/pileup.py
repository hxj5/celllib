
# Note: now it only supports 10x data.

import getopt
import multiprocessing
import os
import pickle
import sys
import time

from .config import APP
from .plp.config import Config, \
    CFG_DEBUG, \
    CFG_CELL_TAG, CFG_UMI_TAG, CFG_UMI_TAG_BC, \
    CFG_NPROC,     \
    CFG_MIN_COUNT, CFG_MIN_MAF, \
    CFG_INCL_FLAG, CFG_EXCL_FLAG_UMI, CFG_EXCL_FLAG_XUMI, \
    CFG_MIN_LEN, CFG_MIN_MAPQ
from .plp.core import sp_count
from .plp.thread import ThreadData
from .plp.utils import load_region_from_txt, load_snp_from_vcf, \
    load_snp_from_tsv, merge_mtx, merge_tsv, rewrite_mtx
from .plp.zfile import zopen, ZF_F_GZIP, ZF_F_PLAIN

def prepare_config(conf):
    """Prepare configures for downstream analysis
    @param conf  A Config object.
    @return      0 if success, -1 otherwise.
    @note        This function should be call after cmdline is parsed.
    """
    func = "prepare_config"
    if conf.sam_fn:
        if not os.path.isfile(conf.sam_fn):
            sys.stderr.write("[E::%s] sam file does not exist '%s'.\n" % (func, conf.sam_fn))
            return(-1)
    else:
        sys.stderr.write("[E::%s] sam file(s) needed!\n" % (func,))
        return(-1)

    if conf.barcode_fn:
        if os.path.isfile(conf.barcode_fn):
            with zopen(conf.barcode_fn, "rt") as fp:
                conf.barcodes = sorted([x.strip() for x in fp])
            if len(set(conf.barcodes)) != len(conf.barcodes):
                sys.stderr.write("[E::%s] duplicate barcodes!\n" % (func, ))
                return(-1)
        else:
            sys.stderr.write("[E::%s] barcode file does not exist '%s'.\n" % (func, conf.barcode_fn))
            return(-1)
    else:       
        conf.barcodes = None
        sys.stderr.write("[E::%s] barcode file needed!\n" % (func,))
        return(-1)

    if not conf.out_dir:
        sys.stderr.write("[E::%s] out dir needed!\n" % func)
        return(-1)
    if not os.path.isdir(conf.out_dir):
        os.mkdir(conf.out_dir)
    conf.out_region_fn = os.path.join(conf.out_dir, conf.out_prefix + "region.tsv")
    conf.out_sample_fn = os.path.join(conf.out_dir, conf.out_prefix + "samples.tsv")
    conf.out_ad_fn = os.path.join(conf.out_dir, conf.out_prefix + "AD.mtx")
    conf.out_dp_fn = os.path.join(conf.out_dir, conf.out_prefix + "DP.mtx")
    conf.out_oth_fn = os.path.join(conf.out_dir, conf.out_prefix + "OTH.mtx")

    if conf.region_fn:
        if os.path.isfile(conf.region_fn): 
            conf.reg_list = load_region_from_txt(conf.region_fn, verbose = True)
            if not conf.reg_list:
                sys.stderr.write("[E::%s] failed to load region file.\n" % func)
                return(-1)
            sys.stdout.write("[I::%s] count %d regions in %d single cells.\n" % (func, 
                len(conf.reg_list), len(conf.barcodes)))
        else:
            sys.stderr.write("[E::%s] region file does not exist '%s'.\n" % (func, conf.region_fn))
            return(-1)
    else:
        sys.stderr.write("[E::%s] region file needed!\n" % (func,))
        return(-1)

    if conf.snp_fn:
        if os.path.isfile(conf.snp_fn):
            if conf.snp_fn.endswith(".vcf") or conf.snp_fn.endswith(".vcf.gz") \
                    or conf.snp_fn.endswith(".vcf.bgz"):
                conf.snp_set = load_snp_from_vcf(conf.snp_fn, verbose = True)
            else:
                conf.snp_set = load_snp_from_tsv(conf.snp_fn, verbose = True)
            if not conf.snp_set or conf.snp_set.get_n() <= 0:
                sys.stderr.write("[E::%s] failed to load snp file.\n" % func)
                return(-1)
            else:
                sys.stdout.write("[I::%s] %d SNPs loaded.\n" % (func,
                    conf.snp_set.get_n()))           
        else:
            sys.stderr.write("[E::%s] snp file does not exist '%s'.\n" % (func, conf.snp_fn))
            return(-1)      
    else:
        sys.stderr.write("[E::%s] SNP file needed!\n" % (func,))
        return(-1)

    if conf.cell_tag and conf.cell_tag.upper() == "NONE":
        conf.cell_tag = None
    if conf.cell_tag and conf.barcodes:
        pass       
    elif (not conf.cell_tag) ^ (not conf.barcodes):
        sys.stderr.write("[E::%s] should not specify cell_tag or barcodes alone.\n" % (func, ))
        return(-1)
    else:
        sys.stderr.write("[E::%s] should specify cell_tag and barcodes.\n" % (func, ))
        return(-1)        

    if conf.umi_tag:
        if conf.umi_tag.upper() == "AUTO":
            if conf.barcodes is None:
                conf.umi_tag = None
            else:
                conf.umi_tag = CFG_UMI_TAG_BC
        elif conf.umi_tag.upper() == "NONE":
            conf.umi_tag = None
    else:
        sys.stderr.write("[E::%s] umi tag needed!\n" % (func, ))
        return(-1)

    if conf.barcodes:
        with open(conf.out_sample_fn, "w") as fp:
            fp.write("".join([b + "\n" for b in conf.barcodes]))

    if conf.excl_flag < 0:
        if conf.use_umi():
            conf.excl_flag = CFG_EXCL_FLAG_UMI
        else:
            conf.excl_flag = CFG_EXCL_FLAG_XUMI

    return(0)

def usage(fp = sys.stderr):
    s =  "\n" 
    s += "Usage: %s %s <options>\n" % (APP, COMMAND)  
    s += "\n" 
    s += "Options:\n"
    s += "  -s, --sam STR          Indexed sam/bam/cram file.\n"
    s += "  -O, --outdir DIR       Output directory for sparse matrices.\n"
    s += "  -R, --region FILE      A TSV file listing target regions. The first 4 columns shoud be:\n"
    s += "                         chrom, start, end (both 1-based and inclusive), name.\n"
    s += "  -P, --phasedSNP FILE   A TSV or VCF file listing phased SNPs (i.e., containing phased GT).\n"
    s += "  -b, --barcode FILE     A plain file listing all effective cell barcode.\n"
    s += "  -h, --help             Print this message and exit.\n"
    s += "  -D, --debug INT        Used by developer for debugging [%d]\n" % CFG_DEBUG
    s += "\n"
    s += "Optional arguments:\n"
    s += "  -p, --nproc INT        Number of processes [%d]\n" % CFG_NPROC
    s += "  --cellTAG STR          Tag for cell barcodes [%s]\n" % CFG_CELL_TAG
    s += "  --UMItag STR           Tag for UMI, set to None when reads only [%s]\n" % CFG_UMI_TAG
    s += "  --minCOUNT INT         Mininum aggragated count for SNP [%d]\n" % CFG_MIN_COUNT
    s += "  --minMAF FLOAT         Mininum minor allele fraction for SNP [%f]\n" % CFG_MIN_MAF
    s += "  --outputAllReg         If set, output all inputted regions.\n"
    s += "  --countDupHap          If set, UMIs aligned to both haplotypes will be counted\n"
    s += "\n"
    s += "Read filtering:\n"
    s += "  --inclFLAG INT    Required flags: skip reads with all mask bits unset [%d]\n" % CFG_INCL_FLAG
    s += "  --exclFLAG INT    Filter flags: skip reads with any mask bits set [%d\n" % CFG_EXCL_FLAG_UMI
    s += "                    (when use UMI) or %d (otherwise)]\n" % CFG_EXCL_FLAG_XUMI
    s += "  --minLEN INT      Minimum mapped length for read filtering [%d]\n" % CFG_MIN_LEN
    s += "  --minMAPQ INT     Minimum MAPQ for read filtering [%d]\n" % CFG_MIN_MAPQ
    s += "  --countORPHAN     If use, do not skip anomalous read pairs.\n"
    s += "\n"

    fp.write(s)

def show_progress(rv = None):
    return(rv)

def pileup(argv):
    """Core part
    @param argv   A list of cmdline parameters [list]
    @return       0 if success, -1 otherwise [int]
    """
    func = "main"
    ret = -1

    if len(argv) <= 2:
        usage(sys.stderr)
        sys.exit(1)

    start_time = time.time()

    conf = Config()
    opts, args = getopt.getopt(argv[2:], "-s:-O:-R:-P:-b:-h-D:-p:", [
                     "sam=", 
                     "outdir=", 
                     "region=", "phasedSNP=" "barcode=",
                     "help", "debug=",
                     "nproc=", 
                     "cellTAG=", "UMItag=", 
                     "minCOUNT=", "minMAF=", "outputAllReg", "countDupHap",
                     "inclFLAG=", "exclFLAG=", "minLEN=", "minMAPQ=", "countORPHAN"
                ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in   ("-s", "--sam"): conf.sam_fn = val
        elif op in ("-O", "--outdir"): conf.out_dir = val
        elif op in ("-R", "--region"): conf.region_fn = val
        elif op in ("-P", "--phasedsnp"): conf.snp_fn = val
        elif op in ("-b", "--barcode"): conf.barcode_fn = val
        elif op in ("-h", "--help"): usage(sys.stderr); sys.exit(1)
        elif op in ("-D", "--debug"): conf.debug = int(val)

        elif op in ("-p", "--proc"): conf.nproc = int(val)
        elif op in ("--celltag"): conf.cell_tag = val
        elif op in ("--umitag"): conf.umi_tag = val
        elif op in ("--mincount"): conf.min_count = int(val)
        elif op in ("--minmaf"): conf.min_maf = float(val)
        elif op in ("--outputallreg"): conf.output_all_reg = True
        elif op in ("--countduphap"): conf.no_dup_hap = False

        elif op in ("--inclflag"): conf.incl_flag = int(val)
        elif op in ("--exclflag"): conf.excl_flag = int(val)
        elif op in ("--minlen"): conf.min_len = int(val)
        elif op in ("--minmapq"): conf.min_mapq = float(val)
        elif op in ("--countorphan"): conf.no_orphan = False

        else:
            sys.stderr.write("[E::%s] invalid option: '%s'.\n" % (func, op))
            return(-1)

    try:
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
        sys.stdout.write("[I::%s] start time: %s.\n" % (func, time_str))

        cmdline = " ".join(argv)
        sys.stdout.write("[I::%s] CMD: %s\n" % (func, cmdline))

        if prepare_config(conf) < 0:
            raise ValueError("[%s] errcode %d" % (func, -2))
        sys.stderr.write("[I::%s] program configuration:\n" % func)
        conf.show(fp = sys.stderr, prefix = "\t")

        conf.umi_dir = os.path.join(conf.out_dir, "umi")
        if not os.path.exists(conf.umi_dir):
            os.mkdir(conf.umi_dir)

        # extract SNPs for each region
        if conf.debug > 0:
            sys.stderr.write("[D::%s] extract SNPs for each region.\n" % (func,))
        reg_list = []
        for reg in conf.reg_list:
            snp_list = conf.snp_set.fetch(reg.chrom, reg.start, reg.end)
            if snp_list and len(snp_list) > 0:
                reg.snp_list = snp_list
                reg_list.append(reg)
            else:
                if conf.debug > 0:
                    sys.stderr.write("[D::%s] no SNP fetched for region '%s'.\n" % 
                        (func, reg.name))
        sys.stdout.write("[I::%s] %d regions extracted with SNPs.\n" % 
            (func, len(reg_list)))

        if not conf.output_all_reg:
            conf.reg_list = reg_list

        # split region list and save to file       
        m_reg = len(conf.reg_list)
        m_thread = conf.nproc if m_reg >= conf.nproc else m_reg

        reg_fn_list = []
        if m_thread > 1:
            n_reg = m_reg // m_thread
            r_reg = m_reg - n_reg * m_thread
            k_reg = 0
            i_thread = 0
            while k_reg <= m_reg - 1:
                t_reg = n_reg + 1 if i_thread < r_reg else n_reg
                reg_fn = conf.out_prefix + "region.pickle." + str(i_thread)
                reg_fn = os.path.join(conf.out_dir, reg_fn)
                reg_fn_list.append(reg_fn)
                with open(reg_fn, "wb") as fp:
                    pickle.dump(conf.reg_list[k_reg:(k_reg + t_reg)], fp)
                k_reg += t_reg
                i_thread += 1
            for reg in conf.reg_list:  # save memory
                del reg
            conf.reg_list.clear()
            conf.reg_list = None
            conf.snp_set.destroy()
            conf.snp_set = None

        thdata_list = []
        ret_sp = -1
        if m_thread <= 1:
            thdata = ThreadData(
                idx = 0, conf = conf,
                reg_obj = conf.reg_list, is_reg_pickle = False,
                out_region_fn = conf.out_region_fn + ".0",
                out_ad_fn = conf.out_ad_fn + ".0",
                out_dp_fn = conf.out_dp_fn + ".0",
                out_oth_fn = conf.out_oth_fn + ".0",
                out_fn = None
            )
            thdata_list.append(thdata)
            if conf.debug > 0:
                sys.stderr.write("[D::%s] data of thread-%d before sp_count:\n" % (func, 0))
                thdata.show(fp = sys.stderr, prefix = "\t")
            ret_sp, thdata = sp_count(thdata)
        else:
            pool = multiprocessing.Pool(processes = m_thread)
            mp_result = []
            for i in range(m_thread):
                thdata = ThreadData(
                    idx = i, conf = conf,
                    reg_obj = reg_fn_list[i], is_reg_pickle = True,
                    out_region_fn = conf.out_region_fn + "." + str(i),
                    out_ad_fn = conf.out_ad_fn + "." + str(i),
                    out_dp_fn = conf.out_dp_fn + "." + str(i),
                    out_oth_fn = conf.out_oth_fn + "." + str(i),
                    out_fn = None
                )
                thdata_list.append(thdata)
                if conf.debug > 0:
                    sys.stderr.write("[D::%s] data of thread-%d before sp_count:\n" % (func, i))
                    thdata.show(fp = sys.stderr, prefix = "\t")
                mp_result.append(pool.apply_async(
                        func = sp_count, 
                        args = (thdata, ), 
                        callback = show_progress))   # TODO: error_callback?
            pool.close()
            pool.join()
            mp_result = [res.get() for res in mp_result]
            retcode_list = [item[0] for item in mp_result]
            thdata_list = [item[1] for item in mp_result]
            if conf.debug > 0:
                sys.stderr.write("[D::%s] returned values of multi-processing:\n" % func)
                sys.stderr.write("\t%s\n" % str(retcode_list))

        # check running status of each sub-process
        for thdata in thdata_list:         
            if conf.debug > 0:
                sys.stderr.write("[D::%s] data of thread-%d after sp_count:\n" % (func, thdata.idx))
                thdata.show(fp = sys.stderr, prefix = "\t")
            if thdata.ret < 0:
                raise ValueError("[%s] errcode %d" % (func, -3))

        if m_thread <= 1:
            if ret_sp < 0:
                raise ValueError("[%s] errcode %d" % (func, -5))

            thdata = thdata_list[0]
            os.rename(thdata.out_region_fn, conf.out_region_fn)

            if rewrite_mtx(thdata.out_ad_fn, ZF_F_GZIP, 
                           conf.out_ad_fn, "wb", ZF_F_PLAIN, 
                           thdata.nr_reg, len(conf.barcodes), thdata.nr_ad,
                           remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -7))

            if rewrite_mtx(thdata.out_dp_fn, ZF_F_GZIP, 
                           conf.out_dp_fn, "wb", ZF_F_PLAIN, 
                           thdata.nr_reg, len(conf.barcodes), thdata.nr_dp,
                           remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -9)) 

            if rewrite_mtx(thdata.out_oth_fn, ZF_F_GZIP, 
                           conf.out_oth_fn, "wb", ZF_F_PLAIN, 
                           thdata.nr_reg, len(conf.barcodes), thdata.nr_oth,
                           remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -11))
        else:
            if merge_tsv([td.out_region_fn for td in thdata_list], ZF_F_GZIP, 
                         conf.out_region_fn, "wb", ZF_F_PLAIN, 
                         remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -15))

            nr_reg_list = [td.nr_reg for td in thdata_list]

            if merge_mtx([td.out_ad_fn for td in thdata_list], ZF_F_GZIP, 
                         conf.out_ad_fn, "w", ZF_F_PLAIN,
                         nr_reg_list, len(conf.barcodes), sum([td.nr_ad for td in thdata_list]),
                         remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -17))

            if merge_mtx([td.out_dp_fn for td in thdata_list], ZF_F_GZIP, 
                         conf.out_dp_fn, "w", ZF_F_PLAIN,
                         nr_reg_list, len(conf.barcodes), sum([td.nr_dp for td in thdata_list]),
                         remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -19))

            if merge_mtx([td.out_oth_fn for td in thdata_list], ZF_F_GZIP, 
                         conf.out_oth_fn, "w", ZF_F_PLAIN,
                         nr_reg_list, len(conf.barcodes), sum([td.nr_oth for td in thdata_list]),
                         remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -21))

    except ValueError as e:
        sys.stderr.write("[E::%s] '%s'\n" % (func, str(e)))
        sys.stdout.write("[E::%s] Running program failed.\n" % func)
        sys.stdout.write("[E::%s] Quiting ...\n" % func)
        ret = -1

    else:
        sys.stdout.write("[I::%s] All Done!\n" % func)
        ret = 0

    finally:
        sys.stdout.write("[I::%s] CMD: %s\n" % (func, cmdline))

        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        sys.stdout.write("[I::%s] end time: %s\n" % (func, time_str))
        sys.stdout.write("[I::%s] time spent: %.2fs\n" % (func, end_time - start_time))

    return(ret)

def run():
    pileup(sys.argv)

COMMAND = "pileup"

if __name__ == "__main__":
    run()
