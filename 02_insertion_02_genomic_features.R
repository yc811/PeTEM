# Rscript 02_insertion_02_genomic_features.R
start_time <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)

# Bed files
ins_gene_file <- args[1]
ins_CDS_file <- args[2]
ins_5UTR_file <- args[3]
ins_exon_file <- args[4]
ins_intron_file <- args[5]
ins_3UTR_file <- args[6]
ins_promoter_file <- args[7]
ins_IGR_file <- args[8]
TE_file <- args[9]
gene_file <- args[10]
CDS_file <- args[11]
UTR5_file <- args[12]
exon_file <- args[13]
intron_file <- args[14]
UTR3_file <- args[15]
promoter_file <- args[16]
IGR_file <- args[17]
genome_file <- args[18]


#------------------
# Load data
ins_gene <- read.table(ins_gene_file)
ins_CDS <- read.table(ins_CDS_file)
ins_5UTR <- read.table(ins_5UTR_file)
ins_exon <- read.table(ins_exon_file)
ins_intron <- read.table(ins_intron_file)
ins_3UTR <- read.table(ins_3UTR_file)
ins_promoter <- read.table(ins_promoter_file)
ins_IGR <- read.table(ins_IGR_file)
TEnew <- read.table(TE_file)
TEnew$V7 <- TEnew$V3 - TEnew$V2

TE_ins <- function(df) {
  df2 <- df[!duplicated(df[, c("V1", "V2", "V3")]), ]
  return(sum(df2$V3 - df2$V2) + nrow(df2))
}

TE_insertion <- c(TE_ins(ins_gene), TE_ins(ins_CDS), TE_ins(ins_5UTR), TE_ins(ins_exon), 
                  TE_ins(ins_intron), TE_ins(ins_3UTR), TE_ins(ins_promoter), TE_ins(ins_IGR))
TE_insertion2 <- TE_insertion*100/sum(as.numeric(TEnew$V7))

# Load genomic regions
gene <- read.table(gene_file)
CDS <- read.table(CDS_file)
UTR5 <- read.table(UTR5_file)
exon <- read.table(exon_file)
intron <- read.table(intron_file)
UTR3 <- read.table(UTR3_file)
promoter <- read.table(promoter_file)
IGR <- read.table(IGR_file)
genome <- read.table(genome_file)

genome_percent <- function(df){
  df$V7 <- (df$V3 - df$V2) + 1
  sum(df$V7)*100/sum(genome$V3)
}

whole_genome <- c(genome_percent(gene), genome_percent(CDS), genome_percent(UTR5),
                  genome_percent(exon), genome_percent(intron), genome_percent(UTR3),
                  genome_percent(promoter), genome_percent(IGR))
TEinsert_enrich <- log2(TE_insertion2 / whole_genome)

regions <- c("Gene","CDS","5'UTR","Exon","Intron","3'UTR","Promoter","IGR")
colors <- c("Gene"="#000000","CDS"="#333333","5'UTR"="#666666","Exon"="#999999",
            "Intron"="#BBBBBB","3'UTR"="#DDDDDD","Promoter"="#BD7634","IGR"="#D8C5AF")

library(ggplot2)
library(ggbreak)

plot_bar <- function(values, ylab, title, yuplim){
  df <- data.frame(Region=factor(regions, levels=regions), Value=values)
  df$label_y <- ifelse(df$Value < 0, 0, df$Value)  # text y position
  df$vjust <- ifelse(df$Value < 0, -0.5, -0.5)     # adjust text position

  p <- ggplot(df, aes(x=Region, y=Value, fill=Region)) +
    geom_bar(stat="identity", color="black") +
    geom_text(aes(x=Region, y=label_y, label=sprintf("%.2f", Value), vjust=vjust),
              size=6.4, color="firebrick4", fontface="bold") +
    expand_limits(y = max(df$Value) * yuplim) +
    scale_fill_manual(values=colors) +
    theme_classic(base_size=30) +
    theme(axis.text.x = element_text(angle=40, hjust=1),
          axis.ticks.x = element_blank(),
          legend.position="none") +
    labs(x="", y=ylab, title=title)
  return(p)
}

#------------------
png("Inserted_genomic_feature_enrichment.png", width=2000, height=1800, res=300)
plot_bar(TEinsert_enrich, "Log2 enrichment", "TE insertion", 1.8)
dev.off()

png("Inserted_genomic_feature_percentage.png", width=2000, height=1800, res=300)
plot_bar(TE_insertion2, "Percentage (%)", "TE insertion", 1.2)
dev.off()

end_time <- Sys.time()
print(end_time-start_time)


