#!/bin/bash
MainDir=/mnt/sequencedb/PopGen/barbara/simulations/slim_neutral/test/s1

# cd into dir
cd $MainDir

for dir in *; do

sed '/^p/!D' $dir/1_2pops_withOUTGROUP.log | sed '0,/^p1/{//d;}' > $dir/slim_sample.txt

done
