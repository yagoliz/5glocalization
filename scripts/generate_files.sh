#!/bin/bash

python3 ./generate-benchmarks.py -f ../benchmarks/positions/gnb3.json -n 3 -v 10000
python3 ./generate-benchmarks.py -f ../benchmarks/positions/gnb4.json -n 4 -v 10000
python3 ./generate-benchmarks.py -f ../benchmarks/positions/gnb5.json -n 5 -v 10000
python3 ./generate-benchmarks.py -f ../benchmarks/positions/gnb8.json -n 8 -v 10000
python3 ./generate-benchmarks.py -f ../benchmarks/positions/gnb10.json -n 10 -v 10000
python3 ./generate-benchmarks.py -f ../benchmarks/positions/gnb15.json -n 15 -v 10000
python3 ./generate-benchmarks.py -f ../benchmarks/positions/gnb20.json -n 20 -v 10000