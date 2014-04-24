sed '/^p/!D' 1_2pops_withOUTGROUP.log | sed '0,/^p1/{//d;}' > slim_sample.txt
