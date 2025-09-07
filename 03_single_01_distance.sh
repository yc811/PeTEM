#!/bin/bash
set -euo pipefail

LOG="03_single_01_distance.log"
echo "[`date`] Pipeline started" | tee $LOG

start_allall=$(date +%s)

# ====================
# usage
# ====================
usage() {
    echo "Usage: bash 03_single_01_distance.sh -g gene.bed -t TE.bed -eg expression_gene.txt -et expression_TE.txt -lim 15000 -tick 5000 -WD 200" 
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
    -tick) MAJOR_TICK="$2"; shift 2 ;;
    -WD) WD="$2"; shift 2 ;;
    *) usage ;;
  esac
done

echo "[`date`] Input files parsed" | tee -a $LOG

# ====================
# step 1: preprocessing: find out TEs locating near genes
# ====================
echo "[`date`] Step 1. Preprocessing: find out TEs locating near genes" | tee -a $LOG

Rscript - "$GENE_BED" "$TE_BED" "$GENE_EXP" "$TE_EXP"  <<'EOF' 

args <- commandArgs(trailingOnly=TRUE)
gene_bed <- args[1]
TE_bed   <- args[2]
gene_exp <- args[3]
TE_exp   <- args[4]

gene_bed <- read.table(gene_bed, header=F, stringsAsFactors=F)
TE_bed   <- read.table(TE_bed, header=F, stringsAsFactors=F)
gene_exp    <- read.table(gene_exp, header=T, stringsAsFactors=F)
TE_exp    <- read.table(TE_exp, header=T, stringsAsFactors=F)

gene_bed2 <- gene_bed[gene_bed$V4 %in% row.names(gene_exp), ]
TE_bed2   <- TE_bed[TE_bed$V4 %in% row.names(TE_exp), ]

write.table(gene_bed2, "expressed_gene.bed", sep="\t", quote=F, col.names=F, row.names=F)
write.table(TE_bed2, "expressed_TE.bed", sep="\t", quote=F, col.names=F, row.names=F)
EOF

bedtools closest -a expressed_gene.bed -b expressed_TE.bed -id -d -D a > expgene_closest_expTE.bed


# ====================
# step 2: split high/low genes for each stage + generate upstream/downstream beds
# ====================
echo "[`date`] Step 2. Split highly/lowly expressed genes and neighbouring TE regions" | tee -a $LOG

Rscript - "$GENE_BED" "$GENE_EXP" "$LIMIT" <<'EOF' 
args <- commandArgs(trailingOnly=TRUE)

gene_bed <- args[1]
gene_exp <- args[2]
limit <- as.numeric(args[3])

gene_bed <- read.table(gene_bed, header=F, stringsAsFactors=F)
gene_exp    <- read.table(gene_exp, header=T, stringsAsFactors=F)

# function for adjacent regions
adjacent <- function(df, up=TRUE, limit=limit){
  df2 <- df
  if(up){
    # upstream
    df2[df2$V6=="+",3] <- df2[df2$V6=="+",2] - 1
    df2[df2$V6=="+",2] <- df2[df2$V6=="+",3] - (limit-1)
    df2[df2$V6=="-",2] <- df2[df2$V6=="-",3] + 1
    df2[df2$V6=="-",3] <- df2[df2$V6=="-",2] + (limit-1)
  } else {
    # downstream
    df2[df2$V6=="-",3] <- df2[df2$V6=="-",2] - 1
    df2[df2$V6=="-",2] <- df2[df2$V6=="-",3] - (limit-1)
    df2[df2$V6=="+",2] <- df2[df2$V6=="+",3] + 1
    df2[df2$V6=="+",3] <- df2[df2$V6=="+",2] + (limit-1)
  }
  df2$V2[df2$V2 < 1] <- 1
  df2[5] <- 0
  return(df2)
}


# focusing on TEs locating within adjacent regions of genes
clo_TE <- read.table("expgene_closest_expTE.bed", sep="\t", header=F)
clo_TE2 <- clo_TE[(clo_TE$V10)!= ".", ]
clo_TE2 <- clo_TE2[clo_TE2$V13 >= (-limit), ]

gene_exp2 <- gene_exp[row.names(gene_exp) %in% clo_TE2$V4, ]

# highly/lowly expressed genes for each stage
stages <- colnames(gene_exp2)

for(stage in stages){
  vals <- gene_exp2[, stage, drop=FALSE]
  vals <- vals[vals[,1]>0, , drop=FALSE]  # remove unexpressed genes at this stage
  
  # sort and get the 25% highly/lowly expressed genes
  sorted <- vals[order(vals[,1]), , drop=FALSE]
  low  <- head(sorted, nrow(sorted)/4)
  high <- tail(sorted, nrow(sorted)/4)
  
  low_bed  <- merge(gene_bed, low, by.x= "V4", by.y= "row.names", all.y=T)
  high_bed <- merge(gene_bed, high, by.x= "V4", by.y= "row.names", all.y=T)
  low_bed <- low_bed[,c(2,3,4,1,7,6)]
  high_bed <- high_bed[,c(2,3,4,1,7,6)]

  write.table(high_bed, paste0("high_", stage, ".txt"), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(low_bed,  paste0("low_", stage, ".txt"),  row.names=F, col.names=F, quote=F, sep="\t")

  # adjacent up/down
  low_up   <- adjacent(low_bed, up=TRUE,  limit=limit)
  high_up  <- adjacent(high_bed, up=TRUE,  limit=limit)
  low_down <- adjacent(low_bed, up=FALSE, limit=limit)
  high_down<- adjacent(high_bed, up=FALSE, limit=limit)

  write.table(low_up,   paste0("low_",  stage, "_up.bed"),   row.names=F, col.names=F, quote=F, sep="\t")
  write.table(high_up,  paste0("high_", stage, "_up.bed"),   row.names=F, col.names=F, quote=F, sep="\t")
  write.table(low_down, paste0("low_",  stage, "_down.bed"), row.names=F, col.names=F, quote=F, sep="\t")
  write.table(high_down,paste0("high_", stage, "_down.bed"), row.names=F, col.names=F, quote=F, sep="\t")
}
EOF


# ====================
# step 3: process methylation files
# ====================
set -euo pipefail

start_step3=$(date +%s)

echo "[`date`] Step 3. Process methylation files: (1) unzip + filter" | tee -a $LOG

mkdir -p methylation

# Step1: unzip + filter
for f in *.CGmap.gz; do
(
    start=$(date +%s)
    base=${f%.CGmap.gz}
    gunzip -c "$f" | awk '$8>=4 {print $1"\t"$3"\t"$2"\t"$4"\t"$6}' > "methylation/${base}.txt"
    end=$(date +%s)
    echo "[INFO] Preprocessed $f in $((end-start)) sec"
)&
done
wait
echo "[INFO] All replicates done."   | tee -a $LOG

echo "[`date`] Step 3. Process methylation files: (2) calculate average mC of each site at each stage"   | tee -a $LOG

stages=($(ls methylation/*.txt | xargs -n1 basename | cut -d'_' -f1 | sort -u))

for stage in "${stages[@]}"; do
(
    start=$(date +%s)
    echo "[INFO] Processing stage $stage"   | tee -a $LOG
    python3 - <<EOF
import pandas as pd, glob, os, time
stage = "${stage}"
files = glob.glob(f"methylation/{stage}_*.txt")

dfs=[]
for f in files:
    df=pd.read_csv(f,sep="\t",header=None,
                   names=["chr","site","nt","CNN","mC"])
    dfs.append(df)

combined = pd.concat(dfs, axis=0, ignore_index=True)
m = combined.groupby(["chr","site","nt","CNN"], as_index=False)["mC"].mean()
m = m.rename(columns={"mC": f"{stage}_mC"})

m["strand"]=m.nt.replace({"C":"+","G":"-"})
m["name"]=["site_"+str(i+1) for i in range(len(m))]

for t in ["CG","CHG","CHH"]:
    sub=m[m.CNN==t][["chr","site","site","name",f"{stage}_mC","strand"]].dropna()
    sub.to_csv(f"{stage}_{t}.bed",sep="\t",index=False,header=False)
EOF
    end=$(date +%s)
    echo "[INFO] Stage $stage done in $((end-start)) sec"   | tee -a $LOG
)&
done
wait

end_step3=$(date +%s)
echo "[INFO] All stages done."   | tee -a $LOG
echo "[`date`] Step 3. files finished in $((end_step3-start_step3)) sec"   | tee -a $LOG


# ====================
# step 4: intersect with TE / without TE
# ====================
set -euo pipefail

start_step4=$(date +%s)
echo "[`date`] Step 4. TE methylation data" | tee -a $LOG

stages=($(ls *_up.bed | cut -d'_' -f2 | sort -u))  # 所有 stage

for stage in "${stages[@]}"; do
(
    start=$(date +%s)
    echo "[INFO] Processing stage $stage"   | tee -a $LOG

    for expr in low high; do
        for dir in up down; do
            input="${expr}_${stage}_${dir}.bed"

            # 1. with TE
            wTE="wTE_${expr}_${stage}_${dir}.bed"
            bedtools intersect -a "$input" -b expressed_TE.bed > "$wTE"

            for type in CG CHG CHH; do
                bedtools intersect -a "$wTE" -b "${stage}_${type}.bed" -wa -wb > "wTE_${expr}_${stage}_${dir}_${type}.bed"
            done

            # 2. without TE
            #woTE="woTE_${expr}_${stage}_${dir}.bed"
            #bedtools subtract -a "$input" -b expressed_TE.bed > "$woTE"

            #for type in CG CHG CHH; do
            #    bedtools intersect -a "$woTE" -b "${stage}_${type}.bed" -wa -wb > "woTE_${expr}_${stage}_${dir}_${type}.bed"
            #done

            # remove step4 temp files
            rm -f "$wTE" #"$woTE"
        done
    done

    end=$(date +%s)
    echo "[INFO] Stage $stage done in $((end-start)) sec"   | tee -a $LOG
)&
done
wait

end_step4=$(date +%s)
echo "[`date`] Step 4. finished in $((end_step4-start_step4)) sec"   | tee -a $LOG



# ====================
# step 5: plotting
# ====================
set -euo pipefail

start_step5=$(date +%s)
echo "[`date`] Step5: Calculate methylation of the regions with TEs and plot (parallel stages)"   | tee -a $LOG

# all stages
stages=($(ls *_low_*_up*.bed | awk -F'_' '{print $3}' | sort -u))

for stage in "${stages[@]}"; do
{
    start=$(date +%s)
    echo "[INFO] Processing stage $stage"  | tee -a $LOG

    Rscript - "$LIMIT" "$MAJOR_TICK" "$WD" "$stage" <<'EOF'
library(ggplot2)
library(gplots)

args <- commandArgs(trailingOnly=TRUE)
limit <- as.numeric(args[1])
major_tick <- as.numeric(args[2])
WD <- as.numeric(args[3])
stage <- args[4]

#limit <- 15000
#major_tick <- 5000
#WD <- 200      # sliding window size (bp)

# functions
for_up <- function(mydf, mybed){
  mydf2 <- mydf[,c("V4","V8","V11","V12")]
  mybed2 <- mybed[,c("V2","V3","V4","V6")]
  mymy <- merge(mydf2, mybed2, by.x=c("V4","V12"), by.y=c("V4","V6"))
  rev <- mymy[mymy$V12=="-",]
  fow <- mymy[mymy$V12=="+",]
  rev$dist <- rev$V8 - rev$V3
  fow$dist <- fow$V2 - fow$V8
  return(rbind(rev,fow))
}

for_down <- function(mydf, mybed){
  mydf2 <- mydf[,c("V4","V8","V11","V12")]
  mybed2 <- mybed[,c("V2","V3","V4","V6")]
  mymy <- merge(mydf2, mybed2, by.x=c("V4","V12"), by.y=c("V4","V6"))
  rev <- mymy[mymy$V12=="-",]
  fow <- mymy[mymy$V12=="+",]
  rev$dist <- rev$V2 - rev$V8
  fow$dist <- fow$V8 - fow$V3
  return(rbind(rev,fow))
}

linedf <- function(df){
  df2 <- df[,c("V11","dist")]
  df3 <- aggregate(df2$V11, list(df2$dist), mean)
  colnames(df3) <- c("dist","V11")
  df_list <- c()
  for(i in 0:RAN2){
    val <- mean(df3[ (abs(df3$dist)>(i*SS+0.5)) & (abs(df3$dist)<=(i*SS+WD)), ]$V11, na.rm=TRUE)
    df_list <- c(df_list,val)
  }
  return(df_list)
}


# processing
types <- c("CG","CHG","CHH")
exprs <- c("low","high")
dirs <- c("up","down")

RAN <- limit # searching range
WD <- WD     # window size (bp)
SS <- 4  
Start <- WD/2
RAN2 <- (RAN-WD/2)/SS


for(type in types){    # "CG","CHG","CHH"
  for(expr in exprs){  # "low","high"
    for(dir in dirs){  # "up","down"
      mydf <- read.table(paste0("wTE_",expr,"_",stage,"_",dir,"_",type,".bed"), sep="\t")
      mybed <- read.table(paste0(expr,"_",stage,".txt"), sep="\t")
      if(dir=="up"){
        df <- for_up(mydf,mybed)
        df$dist <- 0 - df$dist 
      } else {
        df <- for_down(mydf,mybed)
      }
      df$V11 <- df$V11*100
      assign(paste0(expr,"_",dir),df)
    }
  }

  ### geom_point: merge up/down & low/high ###
  point_low <- rbind(get("low_up"), get("low_down"))
  point_low$expr <- "Lowly expressed genes"
  point_high <- rbind(get("high_up"), get("high_down"))
  point_high$expr <- "Highly expressed genes"
  point_all <- rbind(point_low,point_high)

  ### geom_line: use linedf for smoothing ###
  df_up_low   <- linedf(get("low_up"))
  df_down_low <- linedf(get("low_down"))
  df_up_high  <- linedf(get("high_up"))
  df_down_high<- linedf(get("high_down"))

  up_low_line   <- data.frame(distance=c(seq(-Start, -RAN, -SS)), mC=df_up_low)
  up_high_line  <- data.frame(distance=c(seq(-Start, -RAN, -SS)), mC=df_up_high)
  down_low_line <- data.frame(distance=c(seq(Start, RAN, SS)), mC=df_down_low)
  down_high_line<- data.frame(distance=c(seq(Start, RAN, SS)), mC=df_down_high)

  low_line  <- rbind(up_low_line, down_low_line); low_line$expr <- "Lowly expressed genes"
  high_line <- rbind(up_high_line, down_high_line); high_line$expr <- "Highly expressed genes"
  line_all  <- rbind(low_line,high_line)
  line_all <- line_all[complete.cases(line_all),]

  # check boarder
  check_border <- cbind(low_line, high_line)
  check_up <- check_border[check_border$distance < 0, ]
  check_dn <- check_border[check_border$distance > 0, ]
  border_up <- head(check_up[check_up[,2] < check_up[,5], ], 1)
  border_dn <- head(check_dn[check_dn[,2] < check_dn[,5], ], 1)

  up_title <- if(nrow(border_up)>0) paste0("Upstream border: ", border_up$distance[1], " bp") else ""
  dn_title <- if(nrow(border_dn)>0) paste0("Downstream border: ", border_dn$distance[1], " bp") else ""

  
  # for plot
  gap <- limit/10   # x axis gap between TSS and TTS 

  # upstream
  up_df_point <- point_all[point_all$dist < 0,]
  up_df_line  <- line_all[line_all$distance < 0,]
  up_df_point$dist_shift <- up_df_point$dist
  up_df_line$distance_shift <- up_df_line$distance

  # downstream (shift by gap)
  down_df_point <- point_all[point_all$dist > 0,]
  down_df_line  <- line_all[line_all$distance > 0,]
  down_df_point$dist_shift <- down_df_point$dist + gap
  down_df_line$distance_shift <- down_df_line$distance + gap

  # merge
  df_point_all <- rbind(up_df_point, down_df_point)
  #df_line_all  <- rbind(up_df_line, down_df_line)
  
  breaks <- c(seq(-limit, -100, major_tick),     0,   gap, seq(gap+major_tick, gap+limit, major_tick))
  labels <- c(seq(-limit, -100, major_tick), "TSS", "TTS", seq(major_tick,     limit,     major_tick))

  #plot
  #write.table(low_up,"testtestpng.txt", row.names=F, col.names=F, quote=F, sep="\t")
  png(file=paste0(stage,"_distance_updown_TEm",type,"_HighLowGene.png"), width=5000, height=2000, res=400)

  p<- ggplot() +
  geom_point(df_point_all, mapping=aes(x=dist_shift,y=V11,color=expr), size=0.01, alpha=0.1) +
  geom_line(up_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr)) +
  geom_line(down_df_line, mapping=aes(x=distance_shift,y=mC,color=expr,group=expr)) +
  theme_bw() +
  scale_color_manual(values=c("#B44A53","#509ABC")) +
  scale_x_continuous(breaks=breaks, labels=labels) +
  ggtitle(paste0(up_title, "; ", dn_title)) +
  theme(plot.title = element_text(size=12, face="bold"),
    legend.text=element_text(size=15),
    axis.text=element_text(size=12,face="bold"),
    axis.title=element_text(size=18,face="bold")) +
  labs(x="Distance to gene (bp)", y="TE methylation level (%)")

  print(p)

  dev.off()
}

EOF

    end=$(date +%s)
    echo "[INFO] Stage $stage done in $((end-start)) sec"  | tee -a $LOG
}&
done
wait

end_step5=$(date +%s)
echo "[`date`] Step 5. finished in $((end_step5-start_step5)) sec"  | tee -a $LOG


end_allall=$(date +%s)
echo "[`date`] Pipeline finished in $((end_allall-start_allall)) sec" | tee -a $LOG

