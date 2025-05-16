#!/bin/bash
cd prepare
set -Eeuo pipefail

ref_input="../dream/ref.fasta"
p=50
for shape in "11111111" "1001"
do
    echo "Creating IBF for shape $shape"
    index="s${shape}.ibf"
    mkdir -p "s$shape"
    
    if [ -f "s${shape}/ref.0.header"  ];then   
        rm s${shape}/*
    fi

    valik build "$ref_input" --shape $shape --pattern "$p" --write-out --kmer-count-min 0 --kmer-count-max 254 \
                            --output "$index" --fast --without-parameter-tuning --seg-count 8

    mv ref.*.header s$shape/
    mv ref.*.minimiser s$shape/
done

mkdir -p bins
if [ -f "bins/ref_0.fasta" ];then
    rm bins/*
fi
mv ../dream/ref_*.fasta bins/
