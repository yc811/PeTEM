#!/usr/bin/env Rscript

# Rscript 02_TE_insertion_01_TE_families.R -a TE.bed -i TE_overlap_promoter_u1500_d500.bed -T TE_families.bed 

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


#df1 <- read.table("TE.bed", sep="\t", header=FALSE)
#df2 <- read.table("TE_overlap_promoter_u1500_d500.bed", sep="\t", header=FALSE)

#TE_families <- as.data.frame(read.table("TE_families.bed", sep="\t", header=FALSE))


# Count TE families ----------
df_TE <- data.frame(
  TE  = TE_families$V1,
  df1 = sapply(TE_families$V1, function(cl) sum(grepl(cl, df1$V4))),
  df2 = sapply(TE_families$V1, function(cl) sum(grepl(cl, df2$V4)))
)

print(df_TE)


# Count total number
total_df1 <- sum(df_TE$df1)
total_df2 <- sum(df_TE$df2)

# Count percentage
df_TE <- df_TE %>%
  mutate(
    pc_df1 = df1*100/sum(df1),
    pc_df2 = df2*100/sum(df2)
  )

df_long <- df_TE %>%
  select(TE, pc_df1, pc_df2) %>%
  pivot_longer(cols = starts_with("pc_"), names_to = "type", values_to = "percentage") %>%
  mutate(type = factor(type, levels=c("pc_df1","pc_df2"), labels = c("df1","df2")))

# Calculate inserted family enrichment and p value
df_TE <- df_TE %>%
  mutate(
    enrich = (df2/total_df2) / (df1/total_df1)
  )

p_list <- c()
for(i in 1:nrow(df_TE)){
  # 建立列聯表
  a <- df_TE$df2[i]       # df2 TE count
  b <- df_TE$df1[i]       # df1 TE count
  c <- total_df2 - a
  d <- total_df1 - b
  contingency <- matrix(c(a,c,b,d), nrow=2, byrow=TRUE)
  p_list <- c(p_list, fisher.test(contingency)$p.value)
}

df_TE$pvalue <- format.pval(p_list, digits=3, scientific=TRUE)



# Labels
df_TE <- df_TE %>%
  mutate(
    text = paste0(TE, " (", sprintf("%.2f", enrich), ", p=", pvalue, ")"),
    text_y = 100 - (cumsum(c(0, head(pc_df2, -1))) + pc_df2 / 2)
  )

df_long <- df_TE %>%
  select(TE, pc_df1, pc_df2) %>%
  pivot_longer(cols = starts_with("pc_"), names_to = "type", values_to = "percentage") %>%
  mutate(
    type = factor(type, levels=c("pc_df1","pc_df2"), labels = c("All TEs","Inserted TEs")),
    TE = factor(TE, levels=TE_families$V1) # y 軸順序
  )


# Plot
png(file="Insertion_TE_family_enrichment.png", width=3000, height=2200, res=400)

ggplot(df_long, aes(x = type, y = percentage, alluvium = TE)) +
  geom_alluvium(aes(fill = TE), width = 0.3, alpha = 0.6) +
  geom_stratum(aes(stratum = TE, fill = TE), width = 0.3) +
  geom_text(data = df_TE, aes(label = text, x = 2.3, y = text_y), size = 6, hjust = 0, inherit.aes = FALSE) + 
  scale_x_discrete(limits = c("All TEs", "Inserted TEs", "", "", "")) +
  scale_y_continuous(limits = c(0,100), expand=c(0,0)) +
  scale_fill_hue(c = 50, l = 70, h = c(60, 240)) +
  theme_classic() +
  labs(title = "Enriched families of inserted TEs", x="", y="Percentage", fill="TE Types") +
  theme(
    text = element_text(size=20),
    legend.position = "none",
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

