# Rscript 5_cross_condition_correlation.R --DEG DEG.txt --DETE DETE.txt --unexp "$unexp" 

start_time <- Sys.time()

library(optparse)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(reshape2)
library(stringr)


#---- Option parser ----
option_list = list(
  make_option(c("--DEG"), type="character", help="Gene expression DEG file"),
  make_option(c("--DETE"), type="character", help="TE expression DETE file"),
  make_option(c("--unexp"), type="character", default="n", help="Include unexpressed TEs? (y/n)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

include_unexp <- tolower(opt$unexp) == "y"


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

  # ---- Chi-square test (Q1+Q3 vs Q2+Q4) ----
  group1 <- counts["Q1"] + counts["Q3"]
  group2 <- counts["Q2"] + counts["Q4"]
  chisq_tbl <- matrix(c(group1, group2), nrow=2)
  chisq_res <- chisq.test(chisq_tbl)
  pval_text <- paste0("Chi-squared p = ", signif(chisq_res$p.value, 3))

  caption_text <- paste(eq_text, pval_text, sep="\n")

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
    labs(title=title, caption = caption_text, color="Number of\nneighboring\npoints\n ", x=xlab, y=ylab)
  print(p)

  dev.off()

}

#---- Read expression ----
gene_exp<-read.table(opt$DEG, header=T, sep="\t") 
TE_exp<-read.table(opt$DETE, header=T, sep="\t") 

#---- filter unexpressed ----
gene_stage_cols <- setdiff(colnames(gene_exp), grep("^(logFC_|PValue_|FDR_)", colnames(gene_exp), value=TRUE))
TE_stage_cols   <- setdiff(colnames(TE_exp), grep("^(logFC_|PValue_|FDR_)", colnames(TE_exp), value=TRUE))

gene_exp <- gene_exp[rowSums(gene_exp[, gene_stage_cols, drop=FALSE]) != 0, , drop=FALSE]

if(!include_unexp){
  TE_exp <- TE_exp[rowSums(TE_exp[, TE_stage_cols, drop=FALSE]) != 0, , drop=FALSE]
}


#---- Rename columns ----
colnames(gene_exp) <- gsub("^logFC_", "dexp_", colnames(gene_exp))
colnames(TE_exp)   <- gsub("^logFC_", "dTEexp_", colnames(TE_exp))
colnames(gene_exp) <- gsub("^PValue_", "PV_g_", colnames(gene_exp))
colnames(TE_exp)   <- gsub("^PValue_", "PV_TE_", colnames(TE_exp))
colnames(gene_exp) <- gsub("^FDR_", "FDR_g_", colnames(gene_exp))
colnames(TE_exp)   <- gsub("^FDR_", "FDR_TE_", colnames(TE_exp))

#---- Read methylation ----
CG_TE  <- read.table("Tab_TE_CG.txt", header=TRUE, sep="\t")
CHG_TE <- read.table("Tab_TE_CHG.txt", header=TRUE, sep="\t")
CHH_TE <- read.table("Tab_TE_CHH.txt", header=TRUE, sep="\t")

# follow order of stages in delta expression files
fc_cols <- grep("^dexp_", colnames(gene_exp), value=TRUE)
tmp <- str_match(fc_cols, "^dexp_([A-Za-z0-9]+)_([A-Za-z0-9]+)$")[,2:3]
if(is.null(nrow(tmp))) tmp <- matrix(tmp, nrow=1) 
stages_pairs <- tmp

# calculate delta methylation
for(k in seq_len(nrow(stages_pairs))){
  si <- stages_pairs[k,1]
  sj <- stages_pairs[k,2]
  
  CG_TE[[paste0("dTECG_", si, "_", sj)]]   <- (CG_TE[[si]] - CG_TE[[sj]]) * 100
  CHG_TE[[paste0("dTECHG_", si, "_", sj)]] <- (CHG_TE[[si]] - CHG_TE[[sj]]) * 100
  CHH_TE[[paste0("dTECHH_", si, "_", sj)]] <- (CHH_TE[[si]] - CHH_TE[[sj]]) * 100
}

#---- Merge with promoter overlaps ----
ins_promoter <- read.table("TE_overlap_promoter.bed", header=FALSE)
ins_promoter2 <- ins_promoter[,c("V4","V10")]

ins3 <- merge(gene_exp, ins_promoter2, by.x="row.names", by.y="V10")
colnames(ins3)[1] <- "gene_id"
ins4 <- merge(ins3, TE_exp, by.x="V4", by="row.names")
colnames(ins4)[colnames(ins4)=="V4"] <- "TE_id"

ins_CG  <- merge(ins4, CG_TE,  by.x="TE_id", by.y="ID")
ins_CHG <- merge(ins4, CHG_TE, by.x="TE_id", by.y="ID")
ins_CHH <- merge(ins4, CHH_TE, by.x="TE_id", by.y="ID")

#---- function for Q2/Q4 boxplots ----
plot_gene_TE_box <- function(df_subset, mC_type, si, sj, mode="Q2"){
  if(nrow(df_subset) == 0) return(NULL)
  
  # output gene-TE pair
  out_tbl <- df_subset[,c("gene_id","TE_id")]
  write.table(out_tbl,
              file=paste0("OUTPUT_5_",mode,"_gene_TE_pairs_", mC_type, "_", si, "_", sj,".txt"),
              sep="\t", row.names=FALSE, quote=FALSE)
  
  # expression
  expr_mat <- df_subset[,c("gene_id", paste0(si,".x"), paste0(sj,".x"))]
  expr_mat_long <- melt(expr_mat, id.vars="gene_id", variable.name="stage", value.name="expr")
  expr_mat_long$expr <- log2(expr_mat_long$expr + 0.1)
  expr_mat_long$stage <- gsub("\\.x$", "", expr_mat_long$stage)
  expr_mat_long$stage <- factor(expr_mat_long$stage, levels=c(sj, si)) # stage2, stage1
  
  # methylation
  meth_mat <- df_subset[,c("TE_id", si, sj)]
  meth_mat_long <- melt(meth_mat, id.vars="TE_id", variable.name="stage", value.name="methyl")
  meth_mat_long$methyl <- meth_mat_long$methyl * 100  
  meth_mat_long$stage <- factor(meth_mat_long$stage, levels=c(sj, si))
  
  # scale factor for mapping methylation to expression y-axis
  expr_range <- range(expr_mat_long$expr, na.rm=TRUE)
  meth_range <- range(meth_mat_long$methyl, na.rm=TRUE)
  scale_factor <- diff(expr_range) / diff(meth_range)
  
  # means
  expr_means <- aggregate(expr ~ stage, expr_mat_long, mean)
  meth_means <- aggregate(methyl ~ stage, meth_mat_long, mean)
  
  # plot
  p <- ggplot() +
    geom_boxplot(data=expr_mat_long, aes(x=stage, y=expr),
                 fill="#BFBFBF", width=0.35, position=position_nudge(x=-0.2)) +
    geom_point(data=expr_means, aes(x=stage, y=expr), shape=18, size=6, color="black",
               position=position_nudge(x=-0.2)) +
    geom_text(data=expr_means, aes(x=stage, y=expr, label=sprintf("%.1f", expr)),
              size=6, vjust=-1.2, color="black", position=position_nudge(x=-0.2)) +
    geom_boxplot(data=meth_mat_long, aes(x=stage, y=methyl*scale_factor + expr_range[1] - meth_range[1]*scale_factor),
                 fill="#E5C1AF", width=0.35, position=position_nudge(x=0.2)) +
    geom_point(data=meth_means, aes(x=stage, y=methyl*scale_factor + expr_range[1] - meth_range[1]*scale_factor),
               shape=18, size=6, color="black", position=position_nudge(x=0.2)) +
    geom_text(data=meth_means, aes(x=stage, y=methyl*scale_factor + expr_range[1] - meth_range[1]*scale_factor,
                                   label=sprintf("%.1f", methyl)),
              size=6, vjust=-1.2, color="black", position=position_nudge(x=0.2)) +
    scale_y_continuous(
      name="Expression (log2 RPKM)",
      sec.axis = sec_axis(~ (. - expr_range[1] + meth_range[1]*scale_factor)/scale_factor,
                          name="Methylation level (%)")
    ) +
    scale_x_discrete(name="Stage") +
    theme_bw() +
    theme(axis.title=element_text(face="bold", size=20),
          axis.text=element_text(face="bold", size=18),
          axis.text.x=element_text(face="bold", size=18),
          axis.text.y=element_text(face="bold", size=18))
  
  ggsave(filename=paste0("OUTPUT_5_",mode,"_boxplot_", mC_type, "_", si, "_", sj,".png"), p, width=6, height=5)
  print(p)
}

#---- Scatter plots + Q2/Q4 boxplots ----
for(k in seq_along(fc_cols)){
  x_col <- fc_cols[k]
  si <- stages_pairs[k,1]
  sj <- stages_pairs[k,2]
  pv_col <- paste0("PV_g_", si, "_", sj)

  # gene exp vs TE exp
  df_sub <- subset(ins4,
    abs(ins4[[x_col]]) > 1 & ins4[[pv_col]] < 0.05
  )
  if(nrow(df_sub) > 0){
    plot_delta_scatter(df_sub,
      x_col, paste0("dTEexp_", si, "_", sj),
      paste0("OUTPUT_5_geneexp_TEexp_change_", si, "_", sj, "_scatter.png"),
      paste0("TE and gene expression changes\nbetween ", si, " and ", sj),
      expression(Delta~"Gene expression level (log2 RPKM FC)"),
      expression(Delta~"TE expression level (log2 RPKM FC)")
    )
  }

  # gene exp vs TE mC + Q2/Q4
  for(mC_type in c("CG","CHG","CHH")){
    df_TE <- get(paste0("ins_", mC_type))
    y_col <- paste0("dTE", mC_type, "_", si, "_", sj)
    df_sub <- subset(df_TE,
      abs(df_TE[[x_col]]) > 1 & df_TE[[pv_col]] < 0.05
    )
    if(nrow(df_sub) > 0){
      # scatter plot
      plot_delta_scatter(df_sub,
        x_col, y_col,
        paste0("OUTPUT_5_geneexp_TEm", mC_type, "_change_", si, "_", sj, "_scatter.png"),
        paste0("TE m", mC_type, " and gene expression changes\nbetween ", si, " and ", sj),
        expression(Delta~"Gene expression level (log2 RPKM FC)"),
        expression(Delta~paste("TE methylation level (%)"))
      )
      # Q2 boxplot
      df_Q2 <- subset(df_TE,
        abs(df_TE[[x_col]]) > 1 &
        df_TE[[pv_col]] < 0.05 &
        df_TE[[x_col]] < 0 &
        df_TE[[y_col]] > 0
      )
      plot_gene_TE_box(df_Q2, mC_type, si, sj, mode="Q2")
      # Q4 boxplot
      df_Q4 <- subset(df_TE,
        abs(df_TE[[x_col]]) > 1 &
        df_TE[[pv_col]] < 0.05 &
        df_TE[[x_col]] > 0 &
        df_TE[[y_col]] < 0
      )
      plot_gene_TE_box(df_Q4, mC_type, si, sj, mode="Q4")
    }
  }

  # TE exp vs TE mC
  for(mC_type in c("CG","CHG","CHH")){
    df_TE <- get(paste0("ins_", mC_type))
    x_TE_col <- paste0("dTEexp_", si, "_", sj)
    y_TE_col <- paste0("dTE", mC_type, "_", si, "_", sj)
    pv_TE_col <- paste0("PV_TE_", si, "_", sj)

    df_sub <- subset(df_TE,
      abs(df_TE[[x_TE_col]]) > 1 & df_TE[[pv_TE_col]] < 0.05
    )
    if(nrow(df_sub) > 0){
      plot_delta_scatter(df_sub,
        x_TE_col, y_TE_col,
        paste0("OUTPUT_5_TEexp_TEm", mC_type, "_change_", si, "_", sj, "_scatter.png"),
        paste0("TE m", mC_type, " and TE expression changes\nbetween ", si, " and ", sj),
        expression(Delta~"TE expression level (log2 RPKM FC)"),
        expression(Delta~paste("TE methylation level (%)"))
      )
    }
  }
}

end_time <- Sys.time()
print(end_time-start_time)
