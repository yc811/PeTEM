#!/bin/bash
set -euo pipefail

LOG="03_single_01_distance.log"
echo "[`date`] Pipeline started" > $LOG

# ====================
# usage
# ====================
usage() {
    echo "Usage: $0 -g gene.bed -t TE.bed -eg expression_gene.txt -et expression_TE.txt -lim 15000"
    exit 1
}

# ====================
# parse args
# ====================
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -g) GENE_BED="$2"; shift 2 ;;
    -t) TE_BED="$2"; shift 2 ;;
    -eg) GENE_EXP="$2"; shift 2 ;;
    -et) TE_EXP="$2"; shift 2 ;;
    -lim) LIMIT="$2"; shift 2 ;;
    *) usage ;;
  esac
done

echo "[`date`] Input files parsed" >> $LOG

# ====================
# step 0: preprocessing
# ====================
Rscript - <<EOF
args <- commandArgs(trailingOnly=TRUE)
gene_bed <- "$GENE_BED"
te_bed   <- "$TE_BED"
gene_exp <- "$GENE_EXP"
te_exp   <- "$TE_EXP"

genes <- read.table(gene_bed, header=F, stringsAsFactors=F)
tes   <- read.table(te_bed, header=F, stringsAsFactors=F)
ge    <- read.table(gene_exp, header=T, stringsAsFactors=F)
te    <- read.table(te_exp, header=T, stringsAsFactors=F)

genes_exp <- genes[genes\$V4 %in% ge\$GeneID, ]
tes_exp   <- tes[tes\$V4 %in% te\$TEID, ]

write.table(genes_exp, "expressed_gene.bed", sep="\t", quote=F, col.names=F, row.names=F)
write.table(tes_exp, "expressed_TE.bed", sep="\t", quote=F, col.names=F, row.names=F)
EOF

# ====================
# step 1: split high/low genes for each stage
# ====================
Rscript - <<EOF
gene_exp <- "$GENE_EXP"
limit <- as.numeric("$LIMIT")

exp <- read.table(gene_exp, header=T, stringsAsFactors=F)
stages <- colnames(exp)[-1]

for(stage in stages){
  vals <- exp[[stage]]
  q25 <- quantile(vals, 0.25)
  q75 <- quantile(vals, 0.75)
  high <- exp\$GeneID[vals >= q75]
  low  <- exp\$GeneID[vals <= q25]
  write.table(high, paste0("highGene_", stage, ".txt"), row.names=F, col.names=F, quote=F)
  write.table(low,  paste0("lowGene_", stage, ".txt"),  row.names=F, col.names=F, quote=F)
}
EOF

# ====================
# step 2: generate upstream/downstream beds
# ====================
python3 - <<EOF
import pandas as pd

limit = int("$LIMIT")
genes = pd.read_csv("expressed_gene.bed", sep="\t", header=None)

stages = ["MY","NOD","PR","ST","GI"]
for stage in stages:
    for expType in ["high","low"]:
        ids = pd.read_csv(f"{expType}Gene_{stage}.txt", header=None)[0]
        sel = genes[genes[3].isin(ids)]
        # upstream
        up = sel.copy()
        up[1] = up[1] - limit
        up[2] = sel[1]
        up.to_csv(f"{expType}Gene_{stage}_up.bed", sep="\t", header=False, index=False)
        # downstream
        down = sel.copy()
        down[1] = sel[2]
        down[2] = sel[2] + limit
        down.to_csv(f"{expType}Gene_{stage}_down.bed", sep="\t", header=False, index=False)
EOF

# ====================
# step 3: process methylation files
# ====================
python3 - <<EOF
import sys, gzip, os
import pandas as pd
from collections import defaultdict

indir = "methylation/"
files = [f for f in os.listdir(indir) if f.endswith(".CGmap.gz")]

data = defaultdict(lambda: defaultdict(list))

for f in files:
    with gzip.open(os.path.join(indir,f),'rt') as fin:
        for line in fin:
            chrom, base, pos, ctx, dinuc, ratio, count_mC, count_all = line.strip().split()
            pos = int(pos)
            ratio = float(ratio)
            if int(count_all) < 4: continue
            if ctx not in ["CG","CHG","CHH"]: continue
            data[(chrom,pos)][ctx].append(ratio)

out = {"CG":[],"CHG":[],"CHH":[]}
for (chrom,pos), ctxdict in data.items():
    for ctx in ["CG","CHG","CHH"]:
        if ctx in ctxdict:
            avg = sum(ctxdict[ctx])/len(ctxdict[ctx])
            out[ctx].append([chrom,pos-1,pos,avg])

for ctx in ["CG","CHG","CHH"]:
    pd.DataFrame(out[ctx]).to_csv(f"{ctx}.bed", sep="\t", header=False, index=False)
EOF

# ====================
# step 4: intersect with TE / without TE
# ====================
for ctx in CG CHG CHH; do
  for stage in MY NOD PR ST GI; do
    for exp in high low; do
      for dir in up down; do
        bed="${exp}Gene_${stage}_${dir}.bed"
        out_with="withTE_${stage}_${exp}_${dir}_${ctx}.bed"
        out_wo="woTE_${stage}_${exp}_${dir}_${ctx}.bed"
        bedtools intersect -a $bed -b "$TE_BED" -u > $out_with
        bedtools intersect -a $bed -b "$TE_BED" -v > $out_wo
      done
    done
  done
done

# ====================
# step 5: plotting
# ====================
Rscript - <<EOF
library(ggplot2)

limit <- as.numeric("$LIMIT")
stages <- c("MY","NOD","PR","ST","GI")
contexts <- c("CG","CHG","CHH")

for(stage in stages){
  for(ctx in contexts){
    for(dir in c("up","down")){
      for(te in c("withTE","woTE")){
        for(exp in c("high","low")){
          bedfile <- paste0(te,"_",stage,"_",exp,"_",dir,"_",ctx,".bed")
          if(!file.exists(bedfile)) next
          df <- read.table(bedfile)
          df$dist <- if(dir=="up") df$V2 - df$V2[1] else df$V3 - df$V3[1]
          png(paste0(stage,"_distance_",dir,"_",te,"_",exp,"Gene_",ctx,".png"), width=800,height=600)
          ggplot(df,aes(x=dist,y=V4))+
            geom_point(alpha=0.3,color="dodgerblue")+
            geom_smooth(method="lm",se=FALSE,color="red")+
            labs(title=paste(stage,ctx,te,dir,exp),
                 x=paste0("Distance (",dir,", limit=",limit,")"),
                 y="Methylation level")+
            theme_bw()
          dev.off()
        }
      }
    }
  }
}
EOF

echo "[`date`] Pipeline finished" >> $LOG
