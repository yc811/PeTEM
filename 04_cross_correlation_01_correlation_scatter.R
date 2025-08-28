#!/usr/bin/env Rscript
# Rscript 04_cross_correlation_01_correlation_scatter.R -g expression_gene.txt -t expression_TE.txt 

start_time <- Sys.time()

library(optparse)
library(ggplot2)
library(ggpointdensity)
library(viridis)

#---- Option parser ----
option_list = list(
  make_option(c("-g", "--gene"), type="character", help="Gene expression file"),
  make_option(c("-t", "--TE"), type="character", help="TE expression file")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#---- Helper functions ----
quadrant_counts <- function(df, xx, yy) {
  q1 <- sum(df[[xx]] > 0 & df[[yy]] > 0, na.rm=TRUE)
  q2 <- sum(df[[xx]] < 0 & df[[yy]] > 0, na.rm=TRUE)
  q3 <- sum(df[[xx]] < 0 & df[[yy]] < 0, na.rm=TRUE)
  q4 <- sum(df[[xx]] > 0 & df[[yy]] < 0, na.rm=TRUE)
  c(Q1=q1, Q2=q2, Q3=q3, Q4=q4)
}

lm_eqn_text <- function(df, yy, xx){
  m <- lm(df[[yy]] ~ df[[xx]], df)
  a <- format(coef(m)[1], digits=2)
  b <- format(coef(m)[2], digits=2)
  r2 <- format(summary(m)$r.squared, digits=3)
  eq <- paste0("y = ", a, " + ", b, " * x", ", R² = ", r2)
  return(eq)
}


plot_delta_scatter <- function(df, x, y, fname, title, xlab, ylab){
  counts <- quadrant_counts(df, x, y)
  eq_text <- paste0(lm_eqn_text(df, y, x))
  xlim <- max(abs(df[[x]]), na.rm=TRUE)
  ylim <- max(abs(df[[y]]), na.rm=TRUE)


  png(file=fname, width=3000, height=2500, res=400)
  p <- ggplot(df, aes_string(x=x, y=y)) +
    geom_pointdensity(adjust=0.5, size=2) + 
    scale_color_viridis() +
    geom_hline(yintercept=0, linetype="dashed", color="gray70") +
    geom_vline(xintercept=0, linetype="dashed", color="gray70") +
    geom_smooth(method="lm", se=FALSE, col="black") +
    scale_y_continuous(limits=c(-ylim, ylim)) +
    scale_x_continuous(limits=c(-xlim, xlim)) +
    theme_bw() +
    annotate("text", x=xlim, y=ylim, label=paste("Q1:", counts["Q1"]), size=5, hjust=1, fontface="bold") +
    annotate("text", x=-xlim, y=ylim, label=paste("Q2:", counts["Q2"]), size=5, hjust=0, fontface="bold") +
    annotate("text", x=-xlim, y=-ylim, label=paste("Q3:", counts["Q3"]), size=5, hjust=0, fontface="bold") +
    annotate("text", x=xlim, y=-ylim, label=paste("Q4:", counts["Q4"]), size=5, hjust=1, fontface="bold") +
    theme(legend.text=element_text(size=15),
          legend.title=element_text(size=15, face="bold"),
          axis.text=element_text(size=18, face="bold"),
          axis.title=element_text(size=20, face="bold"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.title=element_text(size=20, face="bold", hjust=0.5),
          plot.caption=element_text(size=10, hjust=0)) +
    labs(title=title, caption = eq_text, color="Number of\nneighboring\npoints\n ", x=xlab, y=ylab)
  print(p)
  dev.off()
}

#---- Read expression ----
gene_exp <- read.table(opt$gene, header=TRUE, sep="\t", row.names=1)
#gene_exp<-read.table("expression_gene.txt", header=T, sep="\t") 
gene_exp <- gene_exp[rowSums(gene_exp) != 0, ]
TE_exp   <- read.table(opt$TE, header=TRUE, sep="\t", row.names=1)
#TE_exp<-read.table("expression_TE.txt", header=T, sep="\t") 
TE_exp   <- TE_exp[rowSums(TE_exp) != 0, ]

stages <- intersect(colnames(gene_exp), colnames(TE_exp))
if(length(stages) < 2) stop("Need at least two common stages")

# delta expression
for(i in 1:(length(stages)-1)){
  for(j in (i+1):length(stages)){
    si <- stages[i]; sj <- stages[j]
    gene_exp[[paste0("dexp_", si, "_", sj)]] <- log2((gene_exp[[si]]+0.1)/(gene_exp[[sj]]+0.1))
    TE_exp[[paste0("dTEexp_", si, "_", sj)]] <- log2((TE_exp[[si]]+0.1)/(TE_exp[[sj]]+0.1))
  }
}

#---- Read methylation ----
CG_TE  <- read.table("TE_CG_tab.txt", header=TRUE, sep="\t")
CHG_TE <- read.table("TE_CHG_tab.txt", header=TRUE, sep="\t")
CHH_TE <- read.table("TE_CHH_tab.txt", header=TRUE, sep="\t")

for(i in 1:(length(stages)-1)){
  for(j in (i+1):length(stages)){
    si <- stages[i]; sj <- stages[j]
    CG_TE[[paste0("dTECG_", si, "_", sj)]]   <- (CG_TE[[si]] - CG_TE[[sj]])*100
    CHG_TE[[paste0("dTECHG_", si, "_", sj)]] <- (CHG_TE[[si]] - CHG_TE[[sj]])*100
    CHH_TE[[paste0("dTECHH_", si, "_", sj)]] <- (CHH_TE[[si]] - CHH_TE[[sj]])*100
  }
}

#---- Merge with promoter overlaps ----
ins_promoter <- read.table("TE_overlap_promoter_u1500_d500.bed", header=FALSE)
ins_promoter2 <- ins_promoter[,c("V4","V10")]

ins3 <- merge(gene_exp, ins_promoter2, by.x="row.names", by.y="V10")
colnames(ins3)[1] <- "gene_id"
ins4 <- merge(ins3, TE_exp, by.x="V4", by="row.names")
colnames(ins4)[colnames(ins4)=="Row.names"] <- "TE_id"

ins_CG  <- merge(ins4, CG_TE,  by.x="V4", by.y="ID")
ins_CHG <- merge(ins4, CHG_TE, by.x="V4", by.y="ID")
ins_CHH <- merge(ins4, CHH_TE, by.x="V4", by.y="ID")

#---- Scatter plots ----
for(i in 1:(length(stages)-1)){
  for(j in (i+1):length(stages)){
    si <- stages[i]; sj <- stages[j]

    # gene exp vs TE exp
    plot_delta_scatter(ins4,
      paste0("dexp_", si, "_", sj), paste0("dTEexp_", si, "_", sj),
      paste0("geneexp_TEexp_change_", si, "_", sj, "_scatter.png"),
      paste0("TE and gene expression changes\nbetween ", si, " and ", sj),
      expression(Delta~"Gene expression level (log2 RPKM FC)"),
      expression(Delta~"TE expression level (log2 RPKM FC)")
    )

    # gene exp vs TE mC
    plot_delta_scatter(ins_CG,
      paste0("dexp_", si, "_", sj), paste0("dTECG_", si, "_", sj),
      paste0("geneexp_TEmCG_change_", si, "_", sj, "_scatter.png"),
      paste0("TE mCG and gene expression changes\nbetween ", si, " and ", sj),
      expression(Delta~"Gene expression level (log2 RPKM FC)"),
      expression(Delta~"TE CG methylation level (%)")
    )

    plot_delta_scatter(ins_CHG,
      paste0("dexp_", si, "_", sj), paste0("dTECHG_", si, "_", sj),
      paste0("geneexp_TEmCHG_change_", si, "_", sj, "_scatter.png"),
      paste0("TE mCHG and gene expression changes\nbetween ", si, " and ", sj),
      expression(Delta~"Gene expression level (log2 RPKM FC)"),
      expression(Delta~"TE CHG methylation level (%)")
    )

    plot_delta_scatter(ins_CHH,
      paste0("dexp_", si, "_", sj), paste0("dTECHH_", si, "_", sj),
      paste0("geneexp_TEmCHH_change_", si, "_", sj, "_scatter.png"),
      paste0("TE mCHH and gene expression changes\nbetween ", si, " and ", sj),
      expression(Delta~"Gene expression level (log2 RPKM FC)"),
      expression(Delta~"TE CHH methylation level (%)")
    )

    # TE exp vs TE mC
    plot_delta_scatter(ins_CG,
      paste0("dTEexp_", si, "_", sj), paste0("dTECG_", si, "_", sj),
      paste0("TEexp_TEmCG_change_", si, "_", sj, "_scatter.png"),
      paste0("TE mCG and TE expression changes\nbetween ", si, " and ", sj),
      expression(Delta~"TE expression level (log2 RPKM FC)"),
      expression(Delta~"TE CG methylation level (%)")
    )

    plot_delta_scatter(ins_CHG,
      paste0("dTEexp_", si, "_", sj), paste0("dTECHG_", si, "_", sj),
      paste0("TEexp_TEmCHG_change_", si, "_", sj, "_scatter.png"),
      paste0("TE mCHG and TE expression changes\nbetween ", si, " and ", sj),
      expression(Delta~"TE expression level (log2 RPKM FC)"),
      expression(Delta~"TE CHG methylation level (%)")
    )

    plot_delta_scatter(ins_CHH,
      paste0("dTEexp_", si, "_", sj), paste0("dTECHH_", si, "_", sj),
      paste0("TEexp_TEmCHH_change_", si, "_", sj, "_scatter.png"),
      paste0("TE mCHH and TE expression changes\nbetween ", si, " and ", sj),
      expression(Delta~"TE expression level (log2 RPKM FC)"),
      expression(Delta~"TE CHH methylation level (%)")
    )
  }
}


end_time <- Sys.time()
print(end_time-start_time)

