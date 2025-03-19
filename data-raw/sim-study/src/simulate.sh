#!/bin/bash

IFS1=$'\n' read -d '' -r -a ivec < par/i.txt
IFS2=$'\n' read -d '' -r -a nvec < par/n.txt
IFS3=$'\n' read -d '' -r -a rvec < par/r.txt
IFS4=$'\n' read -d '' -r -a kvec < par/k.txt

nchains=0
for i in "${ivec[@]}"; do
    for n in "${nvec[@]}"; do
        nchains=$((nchains+n))
    done
done

seed=0
for r in "${rvec[@]}"; do
    for k in "${kvec[@]}"; do
                Rscript --vanilla -e "set.seed(${seed}); nbbp::rnbbp(${nchains}, $r, $k) |> cat(file=\"data/r${r}_k${k}.txt\")"
                seed=$((seed+1))
    done
done
