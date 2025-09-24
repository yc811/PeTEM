#!/usr/bin/env bash
set -euo pipefail

#####################################
# Run pipeline interactively
#####################################

echo "Select steps to run (y/n):"
read -p "0. Preprocessing? (y/n): " run0
read -p "1. TE distribution? (y/n): " run1
read -p "2. Promoter-embedded TE families? (y/n): " run2
read -p "3-1. TE impact distance: preprocessing? (y/n): " run3_1
read -p "3-2. TE impact distance: plot? (y/n): " run3_2
read -p "4. Correlation single condition? (y/n): " run4
read -p "5. Correlation across conditions? (y/n): " run5

#####################################
# Collect input files (ask once for shared ones)
#####################################

# BED files shared
if [[ "$run0" == "y" || "$run1" == "y" || "$run3_2" == "y" ]]; then
  read -p "Gene BED file: " gene_bed
  [[ ! -f "$gene_bed" ]] && { echo "Error: $gene_bed not found"; exit 1; }
fi

if [[ "$run0" == "y" || "$run1" == "y" || "$run2" == "y" || "$run3_2" == "y" ]]; then
  read -p "TE BED file: " te_bed
  [[ ! -f "$te_bed" ]] && { echo "Error: $te_bed not found"; exit 1; }
fi

# Genome fasta index file (shared step 0 & step 1)
if [[ "$run0" == "y" ]]; then
  read -p "Genome fasta index: " faidx
  [[ ! -f "$faidx" ]] && { echo "Error: $faidx not found"; exit 1; }
fi

# DEG + DETE files (run0, run3-2, run4, run5) 
if [[ "$run0" == "y" || "$run3_2" == "y" || "$run4" == "y" || "$run5" == "y" ]]; then
  read -p "DEG file: " deg_file
  [[ ! -f "$deg_file" ]] && { echo "Error: $deg_file not found"; exit 1; }
  read -p "DETE file: " dete_file
  [[ ! -f "$dete_file" ]] && { echo "Error: $dete_file not found"; exit 1; }
fi

# Convert DEG + DETE files to expression.txt (run0, run3-2, run4)
if [[ "$run0" == "y" || "$run3_2" == "y" || "$run4" == "y" ]]; then
  Rscript -e '
    extract_expression_only <- function(infile, outfile) {
      df <- read.table(infile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
      keep <- !grepl("^(logFC_|PValue_|FDR_)", colnames(df))
      df_expr <- df[, keep, drop=FALSE]
      write.table(df_expr, file=outfile, sep="\t", quote=FALSE, col.names=NA)
    }
    extract_expression_only("'"$deg_file"'", "gene_expression.txt")
    extract_expression_only("'"$dete_file"'", "TE_expression.txt")
  '
  gene_exp="gene_expression.txt"
  te_exp="TE_expression.txt"
fi

# Include unexpressed TEs in the analyses (run3-2, run4, run5) 
if [[ "$run3_2" == "y" || "$run4" == "y" || "$run5" == "y" ]]; then
  read -p "Include unexpressed TEs? (y/n, default n): " unexp
  unexp=${unexp:-n}
fi

# Promoter upstream/downstream length
if [[ "$run0" == "y" ]]; then
  read -p "Promoter upstream length from TSS (default 1500): " up
  read -p "Promoter downstream length from TSS (default 500): " dn
  up=${up:-1500}
  dn=${dn:-500}
fi

# Methylation (step 0)
if [[ "$run0" == "y" ]]; then
  read -p "Methylation CGmap.gz files (space separated): " -a meth_files
  for f in "${meth_files[@]}"; do
    [[ ! -f "$f" ]] && { echo "Error: $f not found"; exit 1; }
  done
fi

# TE distribution (step 1)
if [[ "$run1" == "y" ]]; then
  read -p "CDS BED file: " cds_bed
  [[ ! -f "$cds_bed" ]] && { echo "Error: $cds_bed not found"; exit 1; }
  read -p "Exon BED file: " exon_bed
  [[ ! -f "$exon_bed" ]] && { echo "Error: $exon_bed not found"; exit 1; }
  read -p "5'UTR BED file: " utr5_bed
  [[ ! -f "$utr5_bed" ]] && { echo "Error: $utr5_bed not found"; exit 1; }
  read -p "3'UTR BED file: " utr3_bed
  [[ ! -f "$utr3_bed" ]] && { echo "Error: $utr3_bed not found"; exit 1; }
fi

# Promoter-embedded TE families (step 2)
if [[ "$run2" == "y" ]]; then
  read -p "TE family file: " te_family
  [[ ! -f "$te_family" ]] && { echo "Error: $te_family not found"; exit 1; }
fi


# Step 3-2 specific parameters
if [[ "$run3_2" == "y" ]]; then
  read -p "Limit up-/down-stream range (bp)(e.g. 15000): " limit
  read -p "Tick size (bp)(e.g. 5000): " tick
  read -p "Window size (bp)(e.g. 200): " window
  limit=${limit:-15000}
  tick=${tick:-5000}
  window=${window:-200}
fi

# Step 4 specific parameters
if [[ "$run4" == "y" ]]; then
  read -p "Window number (default 156): " wd_num
  read -p "y-axis limit for gene expression vs TE/promoter mC plot (CG, default 50): " ylim_cg
  read -p "y-axis limit for gene expression vs TE/promoter mC plot (CHG, default 10): " ylim_chg
  read -p "y-axis limit for gene expression vs TE/promoter mC plot (CHH, default 10): " ylim_chh
  read -p "y-axis limit for TE expression vs TE mC plot (CH, default 15): " ylim_te_ch
  read -p "y-axis limit for TE expression vs TE mC plot (CG, default 30): " ylim_te_cg
  wd_num=${wd_num:-156}
  ylim_cg=${ylim_cg:-50}
  ylim_chg=${ylim_chg:-10}
  ylim_chh=${ylim_chh:-10}
  ylim_te_ch=${ylim_te_ch:-15}
  ylim_te_cg=${ylim_te_cg:-30}
fi


#####################################
# Execute steps
#####################################

# Step 0 + 3-1 parallel if both selected
if [[ "$run0" == "y" || "$run3_1" == "y" ]]; then
    if [[ "$run0" == "y" && "$run3_1" == "y" ]]; then
        bash 0_preprocessing.sh -g "$gene_bed" -t "$te_bed" -eg "$gene_exp" -et "$te_exp" -fai "$faidx" -up "$up" -dn "$dn" -m "${meth_files[@]}" &
        pid0=$!
        bash 3_1_TE_impact_distance_preprocess.sh &
        pid3_1=$!
        wait $pid0
        wait $pid3_1
    else
        [[ "$run0" == "y" ]] && bash 0_preprocessing.sh -g "$gene_bed" -t "$te_bed" -eg "$gene_exp" -et "$te_exp" -fai "$faidx" -up "$up" -dn "$dn" -m "${meth_files[@]}"
        [[ "$run3_1" == "y" ]] && bash ../3_1_TE_impact_distance_preprocess.sh
    fi
fi

# Step 1
[[ "$run1" == "y" ]] && bash 1_TE_distribution.sh \
  -g "$gene_bed" -c "$cds_bed" \
  -5 "$utr5_bed" -e "$exon_bed" -3 "$utr3_bed" \
  -p promoter.bed -t "$te_bed" -fai "$faidx"

# Step 2
[[ "$run2" == "y" ]] && Rscript 2_TE_families.R -a "$te_bed" -i TE_overlap_promoter.bed -T "$te_family"

# Step 3-2
[[ "$run3_2" == "y" ]] && bash 3_2_TE_impact_distance_plot.sh -g "$gene_bed" -t "$te_bed" -eg "$gene_exp" -et "$te_exp" -lim "$limit" -tick "$tick" -WD "$window" -unexp "$unexp" 

# Step 4
[[ "$run4" == "y" ]] && Rscript 4_single_condition_correlation.R --eg "$gene_exp" --et "$te_exp" --unexp "$unexp" --wd_num "$wd_num" --ylim_CG "$ylim_cg" --ylim_CHG "$ylim_chg" --ylim_CHH "$ylim_chh" --ylim_TEexpTEmC_CH "$ylim_te_ch" --ylim_TEexpTEmC_CG "$ylim_te_cg"

# Step 5 (use raw DEG/DETE directly)
[[ "$run5" == "y" ]] && Rscript 5_cross_condition_correlation.R --DEG "$deg_file" --DETE "$dete_file" --unexp "$unexp" 

echo "[DONE] Pipeline completed successfully."
