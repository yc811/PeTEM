#!/bin/bash

set -euo pipefail

LOG="LOG_3_1_TE_impact_distance_preprocess.log"

start_step=$(date +%s)

echo "[`date`] Preprocess methylation files: (1) unzip + filter" | tee -a $LOG

rm -rf pre_step3
mkdir -p pre_step3

# Step1: unzip + filter
for f in *.CGmap.gz; do
(
    start=$(date +%s)
    base=${f%.CGmap.gz}
    gunzip -c "$f" | awk '$8>=4 {print $1"\t"$3"\t"$2"\t"$4"\t"$6}' > "pre_step3/${base}.txt"
    end=$(date +%s)
    echo "[INFO] Preprocessed $f in $((end-start)) sec"
)&
done
wait
echo "[INFO] All replicates done."   | tee -a $LOG

echo "[`date`] Preprocess methylation files: (2) calculate average mC of each site at each stage"   | tee -a $LOG

stages=($(ls pre_step3/*.txt | xargs -n1 basename | cut -d'_' -f1 | sort -u))

for stage in "${stages[@]}"; do
(
    start=$(date +%s)
    echo "[INFO] Processing stage $stage"   | tee -a $LOG
    python3 - <<EOF
import pandas as pd, glob, os, time
stage = "${stage}"
files = glob.glob(f"pre_step3/{stage}_*.txt")

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

for type in ["CG","CHG","CHH"]:
    sub=m[m.CNN==type][["chr","site","site","name",f"{stage}_mC","strand"]].dropna()
    sub.to_csv(f"pre_step3/pre3_{stage}_{type}.bed",sep="\t",index=False,header=False)
EOF
    end=$(date +%s)
    echo "[INFO] Stage $stage done in $((end-start)) sec"   | tee -a $LOG
)&
done
wait

end_step=$(date +%s)
echo "[INFO] All stages done."   | tee -a $LOG
echo "[`date`] Preprocessed methylation files finished in $((end_step-start_step)) sec"   | tee -a $LOG

