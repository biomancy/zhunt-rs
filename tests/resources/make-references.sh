#!/bin/bash
set -euxo pipefail

gcc zhunt-reference.c -Og -o zhunt -lm

./zhunt 12 6 12 HSV-1.NC_001806.2-crop.txt &
./zhunt 8 1 8 mtDNA.NC_012920.1.txt &
./zhunt 8 1 8 test.txt &
./zhunt 11 3 11 chrY.NC_000024.10-crop.txt &
./zhunt 8 3 8 random.txt &
./zhunt 12 3 12 short.txt &
wait

rm zhunt
