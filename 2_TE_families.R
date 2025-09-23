#!/usr/bin/env Rscript

# Rscript 2_TE_families.R -a TE.bed -i TE_overlap_promoter.bed -T TE_family.txt 

start_time <- Sys.time()

# ==== Load packages ====
library(optparse)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggalluvial)


# Options -------------
option_list <- list(
  make_option(c("-a", "--all"),   type="character", help="All TE annotation file (bed format)"),
  make_option(c("-i", "--ins"),   type="character", help="Annotation file of TEs overlapping with promoters (bed format)"),
  make_option(c("-T", "--TE"),    type="character", help="TE families file (bed format)")
)

args <- parse_args(OptionParser(option_list=option_list))

if (is.null(args$all) | is.null(args$ins) | is.null(args$TE)) {
  cat("Missing required arguments\n\n")
  print_help(OptionParser(option_list=option_list))
  quit(status=1)
}


# Read input files ----------
df1 <- read.table(args$all, sep="\t", header=FALSE)
df2 <- read.table(args$ins, sep="\t", header=FALSE)
TE_families <- as.data.frame(read.table(args$TE, sep="\t", header=FALSE))


df1<-df1[,c(1:6)]
df1<-df1[!duplicated(df1),]

df2<-df2[,c(1:6)]
df2<-df2[!duplicated(df2),]

# Count TE families ----------
df1m <- merge(df1, TE_families, by.x="V4", by.y="V1")
df2m <- merge(df2, TE_families, by.x="V4", by.y="V1")

c1 <- as.data.frame(table(df1m$V2.y))
c2 <- as.data.frame(table(df2m$V2.y))

df_TE <- merge(c1, c2, by="Var1", all=TRUE)
df_TE[is.na(df_TE)] <- 0
colnames(df_TE)<-c("TE", "df1", "df2")

print(df_TE)

# totals
total_df1 <- sum(df_TE$df1)
total_df2 <- sum(df_TE$df2)

# percentage
df_TE$pc_df1 <- df_TE$df1 * 100 / total_df1
df_TE$pc_df2 <- df_TE$df2 * 100 / total_df2

# enrichment
df_TE$enrich <- (df_TE$df2 / total_df2) / (df_TE$df1 / total_df1)

# Fisher exact test
p_list <- numeric(nrow(df_TE))
for(i in seq_len(nrow(df_TE))){
  a <- df_TE$df2[i]
  b <- df_TE$df1[i]
  c <- total_df2 - a
  d <- total_df1 - b
  contingency <- matrix(c(a,c,b,d), nrow=2, byrow=TRUE)
  p_list[i] <- fisher.test(contingency)$p.value
}

df_TE$pvalue_num <- p_list     
df_TE$pvalue <- format.pval(p_list, digits=3, scientific=TRUE)

write.table(df_TE, file="Promoter_embedded_TE_family.txt", sep="\t", quote=F, row.names=T)

# labels
df_TE$text <- paste0(df_TE$TE, " (", sprintf("%.2f", df_TE$enrich), ", p=", df_TE$pvalue, ")")
df_TE$text_y <- 100-(cumsum(df_TE$pc_df2) - df_TE$pc_df2/2)

# long format (不用 tidyr)
df_long <- data.frame(
  TE = rep(df_TE$TE, times=2),
  type = rep(c("All TEs","Inserted TEs"), each=nrow(df_TE)),
  percentage = c(df_TE$pc_df1, df_TE$pc_df2),
  stringsAsFactors=FALSE
)


# enrich > 1，get 3 largest enrichment & smallest pvalue 
df_pos <- df_TE[df_TE$enrich > 1, ]
df_pos <- df_pos[order(df_pos$pvalue_num, -df_pos$enrich), ][1:3, ]

# enrich < 1，get 3 smallest enrichment & smallest pvalue 
df_neg <- df_TE[df_TE$enrich < 1, ]
df_neg <- df_neg[order(df_neg$pvalue_num, df_neg$enrich), ][1:3, ]

df_label <- rbind(df_pos, df_neg)


# modify the labels
df_label <- df_label[order(df_label$text_y), ] 
df_label$y_start <- df_label$text_y
df_label$text_y <- seq(25, 75, length.out = nrow(df_label))
df_label$y_end <- df_label$text_y

# Plot
png(file="Promoter_embedded_TE_family_enrichment.png", width=4500, height=2200, res=400)

ggplot(df_long, aes(x = type, y = percentage, alluvium = TE)) +
  geom_alluvium(aes(fill = TE), width = 0.3, alpha = 0.6) +
  geom_stratum(aes(stratum = TE, fill = TE), width = 0.3) +

  geom_segment(data = df_label,
               aes(x = 2.12, xend = 2.5, y = y_start, yend = y_end),
               color = "gray30", linewidth = 0.7, inherit.aes = FALSE) +

  geom_text(data = df_label,
            aes(x = 2.6, y = y_end, label = text),
            hjust = 0, size = 6, inherit.aes = FALSE) +
  scale_x_discrete(limits = c("All TEs", "Inserted TEs", "", "", "")) +
  scale_y_continuous(limits = c(0,100), expand=c(0,0)) +
  scale_fill_hue(c = 50, l = 70, h = c(10, 280)) +
  theme_classic() +
  labs(title = "Enriched families of inserted TEs", x="", y="Percentage", fill="TE Types") +
  theme(
    text = element_text(size=20),
    plot.title = element_text(hjust=0.5),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.text=element_text(size=18, face="bold"),
    axis.line.x.bottom = element_line(color='black'),
    axis.line.y.left = element_line(color='black')
  )

dev.off()

end_time <- Sys.time()
print(end_time-start_time)

