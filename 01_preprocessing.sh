#!/usr/bin/env bash
#####################################
# Input files
#####################################
# Annotation files
### gene.bed --> gene annotation bed file: chr start end name 0 strand
### TE.bed --> TE annotation bed file: chr start end name 0 strand
### TE_families.bed --> one column with all TE families "unique" name
### CDS.bed --> chr start end (name 0 strand)
### 5UTR.bed --> chr start end (name 0 strand)
### exon.bed --> chr start end (name 0 strand)
### 3UTR.bed --> chr start end (name 0 strand)

# Reference genome index
### chrom.size --> chr length

# Methylation data
### stage_replicate.CGmap.gz --> .CGmap.gz files for all the replicates of each stage

#Expression data
### expression_gene.txt --> row name: gene ID; names of each column: stage name; value: average expression across replicates
### expression_TE.txt --> row name: TE ID; names of each column: stage name; value: average expression across replicates


#####################################
set -euo pipefail

#####################################
# Usage
#####################################
usage() {
    echo "Usage: bash 01_preprocessing.sh -g gene.bed -t TE.bed -c chrom.size -m sample1.CGmap.gz [sample2.CGmap.gz ...]"
    exit 1
}

#####################################
# Parse args
#####################################
GENE=""
TE=""
CHROM=""
METH_FILES=()

while [[ $# -gt 0 ]]; do
    case $1 in
        -g) GENE=$2; shift 2;;
        -t) TE=$2; shift 2;;
        -c) CHROM=$2; shift 2;;
        -m)
            shift
            while [[ $# -gt 0 ]] && [[ ! $1 =~ ^- ]]; do
                METH_FILES+=("$1")
                shift
            done
            ;;
        *) usage ;;
    esac
done

[[ -z "$GENE" || -z "$TE" || -z "$CHROM" || ${#METH_FILES[@]} -eq 0 ]] && usage

echo "[INFO] Gene BED: $GENE"
echo "[INFO] TE BED:   $TE"
echo "[INFO] Chrom size: $CHROM"
echo "[INFO] Methylation CGmaps: ${METH_FILES[*]}"

#####################################
# Timer start
#####################################
START=$(date +%s)

#####################################
# 1. Annotations preprocessing
#####################################
echo "[STEP 1] Preprocessing annotations..."

python3 <<'EOF'
import pandas as pd
gene2=pd.read_csv("gene.bed",sep="\t",header=None)
gene3=gene2.copy()
# promoter: -1500 upstream, +500 downstream
gene3.loc[gene3[5]=="+",2]=gene3[1]+499
gene3.loc[gene3[5]=="+",1]=gene3[2]-2000
gene3.loc[gene3[5]=="-",1]=gene3[2]-499
gene3.loc[gene3[5]=="-",2]=gene3[1]+2000
num=gene3._get_numeric_data()
num[num <= 0] = 1
gene3[4]=0
gene3.to_csv("promoter_u1500_d500.bed",sep='\t',index=False,header=False)
EOF

# TE overlap with promoters
bedtools intersect -a "$TE" -b promoter_u1500_d500.bed -wa -wb > TE_overlap_promoter_u1500_d500.bed
# Promoter regions with TE insertions
bedtools intersect -a promoter_u1500_d500.bed -b "$TE" -wa > overlapped_promoter_u1500_d500.bed
bedtools subtract -a overlapped_promoter_u1500_d500.bed -b "$TE" > overlapped_promoterselves_u1500_d500.bed
sort overlapped_promoterselves_u1500_d500.bed | uniq > overlapped_promoterselves_u1500_d500_uniq.bed

# Rename duplicated promoter names
Rscript - <<'EOF'
df <- read.table("overlapped_promoterselves_u1500_d500_uniq.bed", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
df$V4 <- paste0(df$V4, "_", ave(seq_along(df$V4), df$V4, FUN = seq_along))
write.table(df, file="overlapped_promoterselves_u1500_d500_uniq_rename.bed", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
EOF

#####################################
# 2. Methylation preprocessing
#####################################
echo "[STEP 2] Preprocessing methylation data..."

for f in "${METH_FILES[@]}"; do
    echo "[INFO] Processing $f"
    base=$(basename "$f" .CGmap.gz)  # e.g. FB_01

    # CGmap -> wig
    perl CGmap2Wiggle.pl "$f"

    # wig -> bigWig
    #wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig
    #chmod +x wigToBigWig
    ./wigToBigWig ${base}.CG.wig "$CHROM" ${base}.CG.bw
    ./wigToBigWig ${base}.CHG.wig "$CHROM" ${base}.CHG.bw
    ./wigToBigWig ${base}.CHH.wig "$CHROM" ${base}.CHH.bw

    # bigWigAverageOverBed for annotations
    # wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigAverageOverBed
    # chmod +x bigWigAverageOverBed
    for ctx in CG CHG CHH; do
        ./bigWigAverageOverBed ${base}.${ctx}.bw "$TE" ${base}_TE_${ctx}.tab
        ./bigWigAverageOverBed ${base}.${ctx}.bw "$GENE" ${base}_gene_${ctx}.tab
        ./bigWigAverageOverBed ${base}.${ctx}.bw promoter_u1500_d500.bed ${base}_promoter_${ctx}.tab
        ./bigWigAverageOverBed ${base}.${ctx}.bw overlapped_promoterselves_u1500_d500_uniq_rename.bed ${base}_promoterselves_${ctx}.tab
    done
done

#####################################
# 3. Merge in R (per-stage averages)
#####################################
echo "[STEP 3] Merging methylation tables (stage averages) in R..."

Rscript - "$@" <<'EOF'
process_group <- function(feature, ctx){
  tabs <- list.files(pattern=paste0("_",feature,"_",ctx,".tab$"))
  lst <- lapply(tabs, function(f){
    df <- read.table(f, header=F, stringsAsFactors=F)
    # V1 = ID, V3 = depth, V6 = methylation value
    stage_rep <- sub(paste0("_",feature,"_",ctx,".tab$"),"",f)
    parts <- strsplit(stage_rep, "_")[[1]]
    stage <- parts[1]
    replicate <- parts[2]

    data.frame(ID=df$V1,
               value=df$V6,
               depth=df$V3,
               stage=stage,
               replicate = replicate,
               stringsAsFactors=F)
  })

  df <- do.call(rbind, lst)

  # remove features with depth < 5
  df <- df[df$depth >= 5, ]

  # calculate the methylation value of each ID at each stage 
  avg_df <- aggregate(value ~ ID + stage, data=df, FUN=function(x) mean(x, na.rm=TRUE))

  # wide format：ID --> row, stage --> column
  wide <- reshape(avg_df, idvar="ID", timevar="stage", direction="wide")

  # remove "value." in column names
  names(wide) <- sub("^value\\.", "", names(wide))

  # output
  out <- paste0(feature,"_",ctx,"_tab.txt")
  write.table(wide, file=out, sep="\t", quote=F, row.names=F)
  message("Wrote: ", out)
}

for(f in c("TE","gene","promoter","promoterselves")){
  for(ctx in c("CG","CHG","CHH")){
    process_group(f,ctx)
  }
}
EOF


#####################################
# Timer end
#####################################
END=$(date +%s)
RUNTIME=$((END-START))
echo "[DONE] All outputs generated."
echo "[TIME] Total runtime: $RUNTIME seconds ($((RUNTIME/60)) min)"

mkdir wig_bw_tab
mv *.tab *.bw *.wig wig_bw_tab

