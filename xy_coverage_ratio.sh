#!/bin/bash

# modified from https://github.com/edawson/check-sex

xcov=$(echo "scale=4; $(samtools idxstats $1 | grep -P "X\t" | cut -f 3)/$(samtools idxstats $1 | grep -P "X\t" | cut -f 2)" | bc)
ycov=$(echo "scale=4; $(samtools idxstats $1 | grep -P "Y\t" | cut -f 3)/$(samtools idxstats $1 | grep -P "Y\t" | cut -f 2)" | bc)
ratio=$(echo "scale=4; ${xcov}/${ycov}" | bc)

echo
echo "X Coverage: $xcov"
echo "Y Coverage: $ycov"
echo "XY Ratio: $ratio"
echo
echo "Expected Male: 0.9 - 1.5"
echo "Expected Female: 4.5 - 7.0"
echo
