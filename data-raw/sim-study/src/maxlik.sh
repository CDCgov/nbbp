#!/bin/bash

IFS1=$'\n' read -d '' -r -a ivec < par/i.txt
IFS2=$'\n' read -d '' -r -a nvec < par/n.txt
IFS3=$'\n' read -d '' -r -a rvec < par/r.txt
IFS4=$'\n' read -d '' -r -a kvec < par/k.txt

seed=0
for i in "${ivec[@]}"; do
    for n in "${nvec[@]}"; do
        for r in "${rvec[@]}"; do
            for k in "${kvec[@]}"; do
                timeout 10s Rscript --vanilla -e "source(\"src/run_analysis.R\"); fit_ml(\"data/i${i}_n${n}_r${r}_k${k}.txt\", \"maxlik/i${i}_n${n}_r${r}_k${k}.txt\", ${seed})"
                seed=$((seed+1))
            done
        done
    done
done
