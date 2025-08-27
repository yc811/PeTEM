#!/usr/bin/env Rscript
# Rscript 03_single_correlation.R -g expression_gene.txt -t expression_TE.txt --ylim 40

start_time <- Sys.time()

library(zoo)
library(gplots)
library(ggplot2)
library(optparse)

#---- Option parser ----
option_list = list(
  make_option(c("-g", "--gene"), type="character", help="Gene expression file"),
  make_option(c("-t", "--TE"), type="character", help="TE expression file"),
  make_option(c("--ylim"), type="numeric", default=100, help="ylim for gene exp vs TE mC plots")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#---- Functions ----
sort_exp <- function(df, stage_exp){
  df_stage = df[complete.cases(df),]
  stage = df_stage[order(df_stage[,stage_exp]),]
  return(stage)
}

sliding <- function(vec){
  lth <- length(vec)
  step <- floor(lth/100)
  window <- lth-99*step
  return(as.vector(rollapply(zoo(vec*100), width=window, by=step, FUN=mean, align="center")))
}

corr_mC <- function(df, exp, CG, CHG, CHH, method="pearson") {
  data.frame(
    Corr = c(
      cor.test(df[[exp]], df[[CG]], method=method)$estimate,
      cor.test(df[[exp]], df[[CHG]], method=method)$estimate,
      cor.test(df[[exp]], df[[CHH]], method=method)$estimate
    ),
    Pvalue = c(
      cor.test(df[[exp]], df[[CG]], method=method)$p.value,
      cor.test(df[[exp]], df[[CHG]], method=method)$p.value,
      cor.test(df[[exp]], df[[CHH]], method=method)$p.value
    ),
    row.names = c("CG", "CHG", "CHH")
  )
}

corr_exp <- function(df, x, y, method="pearson") {
  corr <- cor.test(df[[x]], df[[y]], method=method)
  data.frame(Corr=corr$estimate, Pvalue=corr$p.value, row.names=paste(x, y, sep="_vs_"))
}

plot_corr_bar <- function(stage, gdf, TEdf, method="pearson") {
  # correlations
  TEexp_TEmC_corr <- corr_mC(TEdf, paste0(stage,"_TEexp"), paste0(stage,"_TE_CG"), paste0(stage,"_TE_CHG"), paste0(stage,"_TE_CHH"), method=method)
  exp_TEmC_corr   <- corr_mC(gdf, paste0(stage,"_exp"), paste0(stage,"_TE_CG"), paste0(stage,"_TE_CHG"), paste0(stage,"_TE_CHH"), method=method)
  exp_TEexp_corr  <- corr_exp(gdf, paste0(stage,"_exp"), paste0(stage,"_TEexp"), method=method)
  
  # combine
  corr_df <- data.frame(
    Comparison = c("TEexp vs TEmCG", "TEexp vs TEmCHG", "TEexp vs TEmCHH",
                   "geneexp vs TEexp", "geneexp vs TEmCG", "geneexp vs TEmCHG", "geneexp vs TEmCHH"),
    Corr = c(
      TEexp_TEmC_corr["CG","Corr"],
      TEexp_TEmC_corr["CHG","Corr"],
      TEexp_TEmC_corr["CHH","Corr"],
      exp_TEexp_corr[1,"Corr"],
      exp_TEmC_corr["CG","Corr"],
      exp_TEmC_corr["CHG","Corr"],
      exp_TEmC_corr["CHH","Corr"]
    ),
    Color = c(rep("#D3A355",3), "#2B5A78", rep("#A4BE78",3))
  )
  corr_df$Comparison <- factor(corr_df$Comparison, levels=rev(corr_df$Comparison))
  
  # plot
  xmin <- min(corr_df$Corr, na.rm=TRUE) - 0.1
  xmax <- max(corr_df$Corr, na.rm=TRUE) + 0.1
  outfile <- paste0(method,"_correlation_bar_",stage,".png")
  png(file=outfile,width=2000,height=1500,res=400)
  print(
    ggplot(corr_df, aes(x=Corr,y=Comparison,fill=Color)) +
      geom_bar(stat="identity") +
      geom_text(aes(label=round(Corr,2)),
                hjust=ifelse(corr_df$Corr>=0, -0.1, 1.1), size=4) +
      scale_fill_identity() +
      xlim(xmin, xmax) +
      labs(x=paste0(method," correlation coefficient"), y=NULL) +
      theme_minimal() +
      theme(
        panel.border=element_rect(colour="black", fill=NA, linewidth=1),
        axis.text=element_text(size=12),
        axis.title.x=element_text(size=14, face="bold")
      )
  )
  dev.off()
}


#---- Read expression ----
gene_exp <- read.table(opt$gene, header=TRUE, sep="\t", row.names=1)
#gene_exp<-read.table("expression_gene.txt", header=T, sep="\t") 
gene_exp <- gene_exp[rowSums(gene_exp) != 0, ]

TE_exp <- read.table(opt$TE, header=TRUE, sep="\t", row.names=1)
#TE_exp<-read.table("expression_TE.txt", header=T, sep="\t") 
TE_exp <- TE_exp[rowSums(TE_exp) != 0, ]

#---- Determine stages ----
stages <- intersect(colnames(gene_exp), colnames(TE_exp))

#---- Read methylation ----
CG_TE <- read.table("TE_CG_tab.txt", header=TRUE, sep="\t")
CHG_TE <- read.table("TE_CHG_tab.txt", header=TRUE, sep="\t")
CHH_TE <- read.table("TE_CHH_tab.txt", header=TRUE, sep="\t")

CG_gene <- read.table("gene_CG_tab.txt", header=TRUE, sep="\t")
CHG_gene <- read.table("gene_CHG_tab.txt", header=TRUE, sep="\t")
CHH_gene <- read.table("gene_CHH_tab.txt", header=TRUE, sep="\t")

CG_promoter <- read.table("promoter_CG_tab.txt", header=TRUE, sep="\t")
CHG_promoter <- read.table("promoter_CHG_tab.txt", header=TRUE, sep="\t")
CHH_promoter <- read.table("promoter_CHH_tab.txt", header=TRUE, sep="\t")

CG_promoterselves <- read.table("promoterselves_CG_tab.txt", header=TRUE, sep="\t")
CHG_promoterselves <- read.table("promoterselves_CHG_tab.txt", header=TRUE, sep="\t")
CHH_promoterselves <- read.table("promoterselves_CHH_tab.txt", header=TRUE, sep="\t")
CG_promoterselves$ID  <- sub("_[0-9]+$", "", CG_promoterselves$ID)
CHG_promoterselves$ID <- sub("_[0-9]+$", "", CHG_promoterselves$ID)
CHH_promoterselves$ID <- sub("_[0-9]+$", "", CHH_promoterselves$ID)


ins_promoter <- read.table("TE_overlap_promoter_u1500_d500.bed", header=FALSE)
ins_promoter2 <- ins_promoter[,c("V4","V10")]

#---- Loop over stages ----
for(stage in stages){

  #---- Merge data ----
  insGene <- merge(ins_promoter2, gene_exp[,stage, drop=FALSE], by.x="V10", by.y="row.names", all.x=TRUE)
  colnames(insGene)[3] <- paste0(stage, "_exp")
  
  insGeneTE <- merge(insGene, TE_exp[,stage, drop=FALSE], by.x="V4", by.y="row.names", all.x=TRUE)
  colnames(insGeneTE)[4] <- paste0(stage, "_TEexp")
  

  # keep only current stage columns for TE methylation
  CG_TE_stage  <- CG_TE[, c("ID", stage), drop=FALSE]
  CHG_TE_stage <- CHG_TE[, c("ID", stage), drop=FALSE]
  CHH_TE_stage <- CHH_TE[, c("ID", stage), drop=FALSE]
  colnames(CG_TE_stage)[2]  <- paste0(stage, "_TE_CG")
  colnames(CHG_TE_stage)[2] <- paste0(stage, "_TE_CHG")
  colnames(CHH_TE_stage)[2] <- paste0(stage, "_TE_CHH")

  insGeneTEmC <- Reduce(function(x,y) merge(x,y,by.x="V4",by.y="ID",all.x=TRUE),
                        list(insGeneTE, CG_TE_stage, CHG_TE_stage, CHH_TE_stage))
  
  # Promoter with TE mC
  CG_promoterselves_stage  <- CG_promoterselves[, c("ID", stage), drop=FALSE]
  CHG_promoterselves_stage <- CHG_promoterselves[, c("ID", stage), drop=FALSE]
  CHH_promoterselves_stage <- CHH_promoterselves[, c("ID", stage), drop=FALSE]
  colnames(CG_promoterselves_stage)[2]  <- paste0(stage, "_promoterselves_CG")
  colnames(CHG_promoterselves_stage)[2] <- paste0(stage, "_promoterselves_CHG")
  colnames(CHH_promoterselves_stage)[2] <- paste0(stage, "_promoterselves_CHH")
  
  insGenePromC <- Reduce(function(x,y) merge(x,y,by.x="V10",by.y="ID",all.x=TRUE),
                         list(insGeneTEmC, CG_promoterselves_stage, CHG_promoterselves_stage, CHH_promoterselves_stage))
 
  # Promoter without TE
  gene_exp_woTE <- gene_exp[!(rownames(gene_exp) %in% ins_promoter$V10), ,drop=FALSE]
  gene_exp_woTE$gene_id <- rownames(gene_exp_woTE)
  
  CG_promoter_stage  <- CG_promoter[, c("ID", stage), drop=FALSE]
  CHG_promoter_stage <- CHG_promoter[, c("ID", stage), drop=FALSE]
  CHH_promoter_stage <- CHH_promoter[, c("ID", stage), drop=FALSE]
  colnames(CG_promoter_stage)[2]  <- paste0(stage, "_promoter_CG")
  colnames(CHG_promoter_stage)[2] <- paste0(stage, "_promoter_CHG")
  colnames(CHH_promoter_stage)[2] <- paste0(stage, "_promoter_CHH")
  
  woTEGenePromC <- Reduce(function(x,y) merge(x,y,by.x="gene_id",by.y="ID",all.x=TRUE),
                          list(gene_exp_woTE, CG_promoter_stage, CHG_promoter_stage, CHH_promoter_stage))
  
  #---- Sorting & sliding ----
  TEdf <- sort_exp(insGeneTEmC[,c("V4","V10",paste0(stage,"_exp"),paste0(stage,"_TEexp"),
                                  paste0(stage,"_TE_CG"),paste0(stage,"_TE_CHG"),paste0(stage,"_TE_CHH"))],
                   paste0(stage,"_TEexp"))
  gdf <- sort_exp(insGeneTEmC[,c("V4","V10",paste0(stage,"_exp"),paste0(stage,"_TEexp"),
                                 paste0(stage,"_TE_CG"),paste0(stage,"_TE_CHG"),paste0(stage,"_TE_CHH"))],
                  paste0(stage,"_exp"))
  
  TE_CG <- sliding(TEdf[[paste0(stage,"_TE_CG")]])
  TE_CHG <- sliding(TEdf[[paste0(stage,"_TE_CHG")]])
  TE_CHH <- sliding(TEdf[[paste0(stage,"_TE_CHH")]])
  
  gTE_CG <- sliding(gdf[[paste0(stage,"_TE_CG")]])
  gTE_CHG <- sliding(gdf[[paste0(stage,"_TE_CHG")]])
  gTE_CHH <- sliding(gdf[[paste0(stage,"_TE_CHH")]])
  gTE_exp <- sliding(gdf[[paste0(stage,"_TEexp")]]/100)
  
  pmC <- sort_exp(insGenePromC[,c("V4","V10",paste0(stage,"_exp"),paste0(stage,"_TEexp"),
                                  paste0(stage,"_promoterselves_CG"),paste0(stage,"_promoterselves_CHG"),paste0(stage,"_promoterselves_CHH"))],
                  paste0(stage,"_exp"))
  pCG <- sliding(pmC[[paste0(stage,"_promoterselves_CG")]])
  pCHG <- sliding(pmC[[paste0(stage,"_promoterselves_CHG")]])
  pCHH <- sliding(pmC[[paste0(stage,"_promoterselves_CHH")]])
  
  wo <- sort_exp(woTEGenePromC[,c("gene_id",paste0(stage),
    paste0(stage,"_promoter_CG"),paste0(stage,"_promoter_CHG"),paste0(stage,"_promoter_CHH"))],
                 paste0(stage))
  woCG <- sliding(wo[[paste0(stage,"_promoter_CG")]])
  woCHG <- sliding(wo[[paste0(stage,"_promoter_CHG")]])
  woCHH <- sliding(wo[[paste0(stage,"_promoter_CHH")]])
  
  #---- Plotting ----
  CG_col <- "#CFA699"
  CHG_col <- "#71A061"
  CHH_col <- "#3871A6"
  TE_col <- "#B9653C"
  pro_wTE_col <- "#7F7F7F"
  pro_woTE_col <- "#BFBFBF"
  
  # TE exp vs TE mC
  png(file=paste0("TEexp_TEmC_line_",stage,".png"), width=2600,height=2200,res=400)
  par(mar=c(5,4.5,4,5)+0.1)
  plot(TE_CHG,lwd=5,lty=1,col=CHG_col,type="l",axes=FALSE,
       ylim=c(0,15),xlim=c(0,100),xlab=NA,ylab=NA,xaxt='n')
  lines(TE_CHH,lwd=5,lty=1,col=CHH_col)
  axis(4,las=1,cex.axis=1.5,font=2)
  mtext(expression(bold("mCH(%)")), side=4, line=3.5, cex=1.5)
  box()
  par(new=TRUE)
  plot(TE_CG,lwd=5,lty=1,col=CG_col,type="l",axes=FALSE,ylim=c(0,30),xlim=c(0,100),xlab=NA,ylab=NA,xaxt='n')
  axis(2, col="black", las=1, cex.axis=1.5, font=2)
  mtext(expression(bold("mCG(%)")), side=2, line=3, cex=1.5)
  mtext("Lowly expressed TEs      Highly expressed TEs", side=1, line=1, cex=1.5, font=2)
  legend("topright", c("CG","CHG","CHH"), text.font=2,bty='n',lty=1,lwd=6,col=c(CG_col,CHG_col,CHH_col),cex=1.8)
  grid(nx=NA, ny=NULL, col="gray70", lty=3, lwd=1)
  dev.off()
  
  # gene exp vs TE exp
  png(file=paste0("geneexp_TEexp_line_",stage,".png"), width=2600,height=2200,res=400)
  par(mar=c(5,4.5,4,5)+0.1)
  plot(gTE_exp,lwd=5,lty=1,col="gray50",type="l",axes=FALSE,xlim=c(0,100),xlab=NA,ylab=NA,xaxt='n')
  axis(2, ylim=c(0,1),col="black",las=1, cex.axis=1.5,font=2)
  mtext(expression(bold("TE expression (log2 RPKM)")),side=2,line=3,cex=1.5)
  mtext("Lowly expressed genes      Highly expressed genes", side=1, line=1, cex=1.5,font=2)
  grid(nx=NA,ny=NULL,col="gray70",lty=3,lwd=1)
  box()
  dev.off()
  
  # gene exp vs TE/promoter mC
  for(mtype in c("CG","CHG","CHH")){
    png(file=paste0("geneexp_TEm",mtype,"_line_",stage,".png"), width=2600,height=2200,res=400)
    par(mar=c(5, 5, 4, 5)+0.1)
    plot(get(paste0("wo",mtype)), lwd=5,lty=1,col=pro_woTE_col,type="l",axes=FALSE,
         ylim=c(0,opt$ylim), xlim=c(0,100), xlab=NA, ylab=NA, xaxt='n')
    lines(get(paste0("p",mtype)), lwd=5,lty=1,col=pro_wTE_col)
    lines(get(paste0("gTE_",mtype)), lwd=5,lty=1,col=CG_col)
    axis(2, las=1, cex.axis=1.5, font=2)
    mtext(expression(bold("Methylation level (%)")), side=2, line=3.5, cex=1.5)
    mtext("Lowly expressed genes      Highly expressed genes", side=1, line=1, cex=1.5,font=2)
    legend("topright", c("TE", "Promoters w TEs", "Promoters w/o TEs"), text.font=2, bty='n', lty=1, lwd=5, col=c(CG_col, pro_wTE_col, pro_woTE_col), cex=1.8)
    box()
    grid(nx=NA, ny=NULL, col="gray70", lty=3, lwd=1)
    dev.off()
  }
  
  # Correlation tables
  # Pearson
  plot_corr_bar(stage, gdf, TEdf, method="pearson")
  # Spearman
  plot_corr_bar(stage, gdf, TEdf, method="spearman")
}

end_time <- Sys.time()
print(end_time-start_time)


