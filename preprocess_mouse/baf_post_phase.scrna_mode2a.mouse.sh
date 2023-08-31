#!/bin/bash
# baf_post_phase.mode2a.mouse.sh 
# - calculate BAF for mouse data, assuming SNPs are called by cellsnp-lite mode 2a,
#   from scRNA-seq data.


function usage() {
    echo
    echo "Usage: $prog [options]"
    echo
    echo "Options:"
    echo "  -N, --name STR       Sample name."
    echo "  -s, --bam FILE       Path to bam file."
    echo "  -b, --barcode FILE   Path to barcode file, one barcode per line."
    echo "  -v, --vcf FILE       Path to phased vcf."
    echo "  -B, --blocks FILE    Path to TSV file containing regions of feature blocks."
    echo "                       If not specified, use built-in feature annotation."
    echo "  -f, --fasta FILE     Path to fasta file."
    echo "  -O, --outdir DIR     Path to output dir."
    echo "  -C, --celltag STR    Cell tag [${def_cell_tag}]"
    echo "  -u, --umi STR        UMI tag. Set to None to count reads [${def_umi}]"
    echo "  -D, --noDUP          If use, duplicate reads will be excluded."
    echo "  -p, --ncores INT     Number of cores [${def_ncores}]"
    echo "  -h, --help           Print this message and exit."
    echo
}


# global settings
work_dir=`cd $(dirname $0); pwd`
prog="baf_post_phase.sh"

ensembl2ucsc=$work_dir/data/ensembl2ucsc.txt
anno_mm10=$work_dir/data/annotate_genes_mm10_cellranger-gex-2020A_20230815.txt


# default settings
def_ncores=1
def_cell_tag=CB
def_umi=UB
use_dup=1
min_count=1
min_maf=0


# check settings
if [ ! -e "$work_dir/utils.sh" ]; then
    echo "Error: utils file $work_dir/utils.sh does not exist!" >&2
    exit 1
fi
source $work_dir/utils.sh

assert_e  "$work_dir/fix_gt.py"  "fix_gt.py"
assert_e  "$ensembl2ucsc"  "ensembl2ucsc file"


# parse args
if [ $# -lt 1 ]; then
    usage
    exit 1
fi

cmdline=`echo $0 $*`

ARGS=`getopt -o N:s:b:v:B:f:O:C:u:Dp:h --long name:,bam:,barcode:,vcf:,blocks:,fasta:,outdir:,celltag:,umi:,noDUP,ncores:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating ..." >&2
    exit 1
fi

eval set -- "$ARGS"
while true; do
    case "$1" in
        -N|--name) sid=$2; shift 2;;
        -s|--bam) bam=$2; shift 2;;
        -b|--barcode) barcode=$2; shift 2;;
        -v|--vcf) vcf=$2; shift 2;;
        -B|--blocks) blocks=$2; shift 2;;
        -f|--fasta) fasta=$2; shift 2;;
        -O|--outdir) out_dir=$2; shift 2;;
        -C|--celltag) cell_tag=$2; shift 2;;
        -u|--umi) umi=$2; shift 2;;
        -D|--noDUP) use_dup=0; shift;;
        -p|--ncores) ncores=$2; shift 2;;
        -h|--help) usage; shift; exit 0;;
        --) shift; break;;
        *) log_err "Internal error!"; exit 1;;
    esac
done

log_msg "CMD: $cmdline"
set -x


# check cmdline args
assert_n  "$sid"  "Sample name"
assert_e  "$bam"  "BAM file"
assert_e  "$barcode"  "Barcode file"
assert_e  "$vcf"  "Phased VCF file"
assert_e  "$fasta"  "FASTA file"

assert_n  "$out_dir"  "Output dir"
if [ ! -e "$out_dir" ]; then
    mkdir -p $out_dir
fi
out_dir=`cd $out_dir; pwd`

if [ -z "$blocks" ]; then
    blocks=$anno_mm10
fi
assert_e  "$blocks"  "Feature file"

if [ -z "$cell_tag" ]; then
    cell_tag=$def_cell_tag
fi

if [ -z "$umi" ]; then
    umi=$def_umi
fi

if [ $use_dup -eq 0 ]; then
    excl_flag=1796
else
    excl_flag=772
fi

if [ -z "$ncores" ]; then
    ncores=$def_ncores
fi


###### Core Part ######

res_dir=$out_dir/result
if [ ! -e "$res_dir" ]; then
    mkdir -p $res_dir
fi

raw_vname=${sid}.vcf.gz
raw_vpath=$vcf

# keep het SNPs only
flt_vname=${raw_vname/.vcf/.het.vcf}
flt_vpath=$res_dir/$flt_vname

log_msg "Keep heterozygous SNPs only."
bcftools view -Oz -i 'GT = "het"' $raw_vpath > $flt_vpath


# add leading chr for chrom names
chr_vname=${flt_vname/.vcf/.chr.vcf}
chr_vpath=$res_dir/$chr_vname

cat $fasta | head -1 | grep -i '^>chr'

if [ $? -eq 0 ]; then     # chrom names have leading chr
    log_msg "Target fasta has leading 'chr'; add 'chr' to VCF chrom names."
    bcftools annotate -Oz --rename-chrs $ensembl2ucsc $lift_vpath > $chr_vpath
else
    log_msg "Target fasta has no leading chr."
    chr_vname=$flt_vname
    chr_vpath=$flt_vpath
fi


# bcftools fixref check
log_msg "(Pre-xcltk-fixref) bcftools fixref check."
bcftools +fixref $chr_vpath -- -f $fasta


# xcltk fixref
fix_vname=${chr_vname/.vcf/.fixref.vcf}
fix_vpath=$res_dir/$fix_vname

log_msg "xcltk fixref."

xcltk fixref -i $chr_vpath  -r $fasta -v |   \
    bgzip -c > $fix_vpath


# bcftools fixref check
log_msg "(Post-xcltk-fixref) bcftools fixref check."
bcftools +fixref $fix_vpath -- -f $fasta


# keep ref-alt heterozygous SNPs only
het_vname=${fix_vname/.vcf/.het.vcf}
het_vpath=$res_dir/$het_vname

log_msg "Keep ref-alt heterozygous SNPs only."
bcftools view -Oz -i 'GT = "0/1" || GT = "1/0"' $fix_vpath > $het_vpath


# fix GT
fgt_vname=${het_vname/.vcf/.fix_gt.vcf}
fgt_vpath=$res_dir/$fgt_vname

log_msg "Fix GT"
python $work_dir/fix_gt.py  -i $het_vpath  -g "0|1"  |   \
    bgzip -c > $fgt_vpath


# filter duplicates and sort
uniq_vname=${fgt_vname/.vcf/.uniq.sort.vcf}
uniq_vpath=$res_dir/$uniq_vname

log_msg "Filter duplicates (chrom + pos) and sort."

zcat $fgt_vpath |    \
    awk '$0 ~ /^#/ {print; next;} ! a[$1":"$2] {print; a[$1":"$2]=1}' |  \
    bcftools sort -Oz > $uniq_vpath


gt_vname=$uniq_vname
gt_vpath=$uniq_vpath


# phase SNPs into haplotype blocks of features
baf_dir=$res_dir/baf

log_msg "Phase SNPs into haplotype blocks of features."

if [ ! -e "$baf_dir" ]; then
    mkdir -p $baf_dir
fi

xcltk pileup -s $bam  -b $barcode  -O $baf_dir  -R $blocks  -P $gt_vpath    \
    -p $ncores                                                 \
    --minCOUNT $min_count  --minMAF $min_maf                   \
    --minLEN 30  --minMAPQ 20  --inclFLAG 0  --exclFLAG $excl_flag   \
    --UMItag $umi  --cellTAG $cell_tag                        \
    --outputAllReg


# move final result
mv $baf_dir $out_dir


###### END ######
log_msg "All Done!"
log_msg "End"

