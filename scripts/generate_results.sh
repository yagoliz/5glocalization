#!/bin/bash

Gs=(3 4 5 10 15 20)
Ws=(1 32 64 128 512)
Ms=(1 2)
for G in ${Gs[@]}; do
    for W in ${Ws[@]}; do
        for M in ${Ms[@]}; do
            ../third_party/tdoapp/build/benchmarks/BenchmarkMean -p 0 -r ../benchmarks/positions/gnb$G.json -w $W -m $M -o ../benchmarks/results/gnb${G}_w${W}_${M}.csv
        done
    done
done
