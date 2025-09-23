# Tutorial
This tutorial explains how to run the PeTEM pipeline using the example files provided in the PeTEM_data folder.

## Installation
Download from github
```bash
git clone https://github.com/yc811/PeTEM.git
cd PeTEM
```
This will download all scripts and the `run_PeTEM.sh` interactive pipeline.

##  Example Data
We provide example input files under [`PeTEM_data/`](PeTEM_data) to test the pipeline.
| File                                                                       | Description                     |
| -------------------------------------------------------------------------- | ------------------------------- |
| `gene.bed`                                                                 | Gene annotation file            |
| `TE.bed`                                                                   | Transposable element annotation |
| `CDS.bed`                                                                  | Coding sequences                |
| `exon.bed`                                                                 | Exon coordinates                |
| `UTR5.bed`                                                                 | 5′ UTR coordinates              |
| `UTR3.bed`                                                                 | 3′ UTR coordinates              |
| `TAIR10.fa.fai`                                                            | Genome index file               |
| `DEG.txt`                                                                  | Differentially expressed genes  |
| `DETE.txt`                                                                 | Differentially expressed TEs    |
| `TE_family.txt`                                                            | TE family annotation            |
| `WT_01.CGmap.gz`, `WT_02.CGmap.gz`, `drdd_01.CGmap.gz`, `drdd_02.CGmap.gz` | Methylation CGmap files         |


## Running the pipeline
PeTEM provides an interactive bash script `run_PeTEM.sh` to select steps and input files.

### Step 1: Launch the interactive pipeline
```
bash run_PeTEM.sh
```
You will be prompted to select which steps to run (enter `y` or `n`):
```
0. Preprocessing? (y/n):
1. TE distribution? (y/n):
2. Promoter-embedded TE families? (y/n):
3-1. TE impact distance: preprocessing? (y/n):
3-2. TE impact distance: plot? (y/n):
4. Correlation single condition? (y/n):
5. Correlation across conditions? (y/n):
```

## Step 2: Input files
After selecting steps, you will be prompted to provide paths to necessary files.
Use the example files in PeTEM_data:
```
Gene BED file: PeTEM_data/gene.bed
TE BED file: PeTEM_data/TE.bed
Genome fasta index: PeTEM_data/TAIR10.fa.fai
DEG file: PeTEM_data/DEG.txt
DETE file: PeTEM_data/DETE.txt
Methylation CGmap.gz files: PeTEM_data/WT_01.CGmap.gz PeTEM_data/WT_02.CGmap.gz PeTEM_data/drdd_01.CGmap.gz PeTEM_data/drdd_02.CGmap.gz
CDS BED file: PeTEM_data/CDS.bed
Exon BED file: PeTEM_data/exon.bed
5' UTR BED file: PeTEM_data/UTR5.bed
3' UTR BED file: PeTEM_data/UTR3.bed
TE family file: PeTEM_data/TE_family.txt
```
For parameters like promoter upstream/downstream length or plotting options, you can press `Enter` to use default values.

## Step 3: Pipeline steps description
| Step | Description                                            |
| ---- | ------------------------------------------------------ |
| 0    | Preprocessing of genome, genes, TEs, methylation files |
| 1    | TE distribution across genomic features                |
| 2    | Promoter-embedded TE families analysis                 |
| 3-1  | TE impact distance preprocessing                       |
| 3-2  | TE impact distance plotting                            |
| 4    | Correlation between TEs and genes (single condition)   |
| 5    | Correlation between TEs and genes across conditions    |

## Step 4: Output files
After successful completion, outputs will include:
* Step 0:
    * Table showing numbers of promoter-embedded TEs and genes with TE embedded in promoters
    * promoter.bed
    * TE_overlap_promoter.bed: Names and location of promoter-embedded TEs and their neighbouring genes
    * Tables showing methylation level of each TE, gene, promoter region under every condition
* Step 1:
    * 01_TE_distribution_enrichment.png
    * 01_TE_distribution_percentage.png
* Step 2:
    * Promoter_embedded_TE_family.txt
    * Promoter_embedded_TE_family_enrichment.png
* Step 3:
    * Figures of TE CG/CHG/CHH methylation around highly/lowly expressed genes under every condition
    * In every figure, it shows the upstream and downstream border (bp) of TE methylation impact, and the number of  highly/lowly expressed genes
* Step 4:
    * Line plots showing correlations between (1) TE methylation vs TE expression, (2) TE methylation, methylation of promoters with TEs/with no TEs vs gene expression, (3) TE expression vs gene expression under each stage. The TE-gene pairs number and the smoothing window size also showed in each figure.
    * Bar plots showing pearson's or spearman's correlation coefficient between (1) TE CG/CHG/CHH methylation vs TE expression, (2) TE CG/CHG/CHH methylation vs gene expression, (3) TE expression vs gene expression under each stage
* Step 5:
    * Scatter plots showing correlations between changes of (1) TE CG/CHG/CHH methylation vs TE expression, (2) TE CG/CHG/CHH methylation vs gene expression, (3) TE expression vs gene expression comparing each two stages. The TE-gene pairs number and the regression line also showed in each figure.
    * Box plots showing the expression and methylation level changes of those negatively correlated gene-TE pairs under each two stages. The name lists of those negatively correlated gene-TE pairs are also in the output.
