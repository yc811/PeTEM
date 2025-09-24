#!/usr/bin/env bash
set -euo pipefail

#####################################
# Usage
#####################################
usage() {
    echo "Usage: bash 0_preprocessing.sh -g gene.bed -t TE.bed -eg expression_gene.txt -et expression_TE.txt -fai genome.fa.fai [-up upstream] [-dn downstream] -m sample1.CGmap.gz [sample2.CGmap.gz ...]"
    echo "  -up   upstream length from TSS (default: 1500)"
    echo "  -dn   downstream length from TSS (default: 500)"
    exit 1
}

#####################################
# Parse args
#####################################
UP=1500  #default
DN=500   #default

GENE=""
TE=""
GENE_EXP=""
TE_EXP=""
faidx=""
METH_FILES=()

while [[ $# -gt 0 ]]; do
    case $1 in
        -g) GENE=$2; shift 2;;
        -t) TE=$2; shift 2;;
        -eg) GENE_EXP=$2; shift 2;;
        -et) TE_EXP=$2; shift 2;;
        -fai) faidx=$2; shift 2;;
        -up) UP=$2; shift 2;;
        -dn) DN=$2; shift 2;;
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

[[ -z "$GENE" || -z "$TE" || -z "$GENE_EXP" || -z "$TE_EXP" || -z "$faidx" || ${#METH_FILES[@]} -eq 0 ]] && usage

echo "[INFO] Gene BED: $GENE"
echo "[INFO] TE BED:   $TE"
echo "[INFO] Gene expression: $GENE_EXP"
echo "[INFO] TE expression:   $TE_EXP"
echo "[INFO] Genome fasta index: $faidx"
echo "[INFO] Methylation CGmaps: ${METH_FILES[*]}"
echo "[INFO] Promoter window: -${UP} upstream, +${DN} downstream"

#####################################
# Timer start
#####################################
START=$(date +%s)

#####################################
# 1. Annotations preprocessing
#####################################

# promoter
echo "[STEP 1] Preprocessing annotations..."

python3 - "$GENE" "$UP" "$DN" <<EOF
import pandas as pd
gene2 = pd.read_csv("${GENE}",sep="\t",header=None)
gene3 = gene2.copy()

up, dn = ${UP}, ${DN}

# promoter: -up upstream, +dn downstream
gene3.loc[gene3[5]=="+", 2] = gene3[1] + dn
gene3.loc[gene3[5]=="+", 1] = gene3[1] - up
gene3.loc[gene3[5]=="-", 1] = gene3[2] - dn
gene3.loc[gene3[5]=="-", 2] = gene3[2] + up

num = gene3._get_numeric_data()
num[num <= 0] = 1
gene3[4] = 0
gene3.to_csv("promoter.bed", sep='\t', index=False, header=False)
EOF

# TE overlap with promoters
bedtools intersect -a "$TE" -b promoter.bed -wa -wb > TE_overlap_promoter.bed

# count TE overlap with promoters
Rscript - "$TE" "$GENE" "$TE_EXP" "$GENE_EXP" <<'EOF'
args <- commandArgs(trailingOnly=TRUE)
te_file <- args[1]
gene_file <- args[2]
te_exp_file <- args[3]
gene_exp_file <- args[4]

te_bed <- read.table(te_file, sep="\t", header=FALSE)
gene_bed <- read.table(gene_file, sep="\t", header=FALSE)
overlap <- read.table("TE_overlap_promoter.bed", sep="\t", header=FALSE)

te_exp <- read.table(te_exp_file, header=TRUE, row.names=1, check.names=FALSE)
gene_exp <- read.table(gene_exp_file, header=TRUE, row.names=1, check.names=FALSE)

te_exp <- te_exp[rowSums(te_exp) != 0, , drop=FALSE]
gene_exp <- gene_exp[rowSums(gene_exp) != 0, , drop=FALSE]

expressed_te <- rownames(te_exp)
expressed_gene <- rownames(gene_exp)

# TE table
all_te <- te_bed$V4
embedded_te <- unique(overlap$V4)

a <- sum(all_te %in% intersect(embedded_te, expressed_te))
b <- sum(all_te %in% setdiff(embedded_te, expressed_te))
c <- sum(all_te %in% setdiff(expressed_te, embedded_te))
d <- sum(all_te %in% setdiff(all_te, union(embedded_te, expressed_te)))

te_tab <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE,
                 dimnames=list(Embedded=c("Yes","No"), Expressed=c("Yes","No")))
te_chi <- chisq.test(te_tab)

# Gene table
all_gene <- gene_bed$V4
embedded_gene <- unique(overlap$V10)

a <- sum(all_gene %in% intersect(embedded_gene, expressed_gene))
b <- sum(all_gene %in% setdiff(embedded_gene, expressed_gene))
c <- sum(all_gene %in% setdiff(expressed_gene, embedded_gene))
d <- sum(all_gene %in% setdiff(all_gene, union(embedded_gene, expressed_gene)))

gene_tab <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE,
                   dimnames=list(Promoter_has_TE=c("Yes","No"), Expressed=c("Yes","No")))
gene_chi <- chisq.test(gene_tab)

# Output txt
sink("OUTPUT_0_embedded_TE_gene_number.txt")
cat("########## TE expression vs promoter embedding ##########\n\n")
print(te_tab)
cat("\nChi-squared p =", te_chi$p.value, "\n\n\n\n")

cat("########## Gene expression vs promoter has TE ##########\n\n")
print(gene_tab)
cat("\nChi-squared p =", gene_chi$p.value, "\n")
sink()
EOF



# Promoter regions with TE insertions
bedtools intersect -a promoter.bed -b "$TE" -wa > overlapped_promoter.bed
bedtools subtract -a overlapped_promoter.bed -b "$TE" > overlapped_promoterselves.bed
sort overlapped_promoterselves.bed | uniq > overlapped_promoterselves_uniq.bed #remove those duplicated regions

# Rename duplicated promoter names
Rscript - <<'EOF'
df <- read.table("overlapped_promoterselves_uniq.bed", sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE)
df$V4 <- paste0(df$V4, "_", ave(seq_along(df$V4), df$V4, FUN = seq_along))
write.table(df, file="overlapped_promoterselves_uniq_rename.bed", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
EOF

#####################################
# 2. Methylation preprocessing
#####################################
echo "[STEP 2] Preprocessing methylation data..."
awk '{print $1"\t"$2}' "$faidx" > chrom.size

for f in "${METH_FILES[@]}"; do
(
    echo "[INFO] Processing $f"
    base=$(basename "$f" .CGmap.gz)  # e.g. FB_01

    # CGmap -> wig
    perl CGmap2Wiggle.pl "$f"

    # wig -> bigWig
    #wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig
    #chmod +x wigToBigWig
    ./wigToBigWig ${base}.CG.wig chrom.size ${base}.CG.bw
    ./wigToBigWig ${base}.CHG.wig chrom.size ${base}.CHG.bw
    ./wigToBigWig ${base}.CHH.wig chrom.size ${base}.CHH.bw

    # bigWigAverageOverBed for annotations
    # wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigAverageOverBed
    # chmod +x bigWigAverageOverBed
    for ctx in CG CHG CHH; do
        ./bigWigAverageOverBed ${base}.${ctx}.bw "$TE" ${base}_TE_${ctx}.tab
        ./bigWigAverageOverBed ${base}.${ctx}.bw "$GENE" ${base}_gene_${ctx}.tab
        ./bigWigAverageOverBed ${base}.${ctx}.bw promoter.bed ${base}_promoter_${ctx}.tab
        ./bigWigAverageOverBed ${base}.${ctx}.bw overlapped_promoterselves_uniq_rename.bed ${base}_promoterselves_${ctx}.tab
    done
) &

done


wait

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
    
    # especially for promoterselves data, merge the methylation values of same promoter
    if(feature == "promoterselves"){
      df$V1 <- sub("_[0-9]+$", "", df$V1)
      df_grouped <- aggregate(. ~ V1, data = df, FUN = sum)
      df_grouped$V6 <- df_grouped$V4 / df_grouped$V3
      df_grouped[is.na(df_grouped)] <- 0
      df <- df_grouped
    }

    data.frame(ID=df$V1,
               value=df$V6,
               depth=df$V3,
               stage=stage,
               replicate = replicate,
               stringsAsFactors=F)
  })

  df <- do.call(rbind, lst)

  # remove features with depth < 5
  C3 <- aggregate(depth ~ ID, data=df, min)
  keep <- C3$ID[C3$depth >= 5]
  df <- df[df$ID %in% keep, ]

  # calculate the methylation value of each ID at each stage 
  avg_df <- aggregate(value ~ ID + stage, data=df, FUN=function(x) mean(x, na.rm=TRUE))

  # wide format：ID --> row, stage --> column
  wide <- reshape(avg_df, idvar="ID", timevar="stage", direction="wide")

  # remove "value." in column names
  names(wide) <- sub("^value\\.", "", names(wide))

  # output
  out <- paste0("Tab_",feature,"_",ctx,".txt")
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

rm overlapped_promoter.bed overlapped_promoterselves.bed overlapped_promoterselves_uniq.bed overlapped_promoterselves_uniq_rename.bed
mkdir wig_bw_tab
mv *.tab *.bw *.wig wig_bw_tab

