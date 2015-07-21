####


## run splitter on slim output, clean file and separate chi9mp, p2 and p3


/mnt/sequencedb/PopGen/barbara/simulations/scripts/perl_test.pl


sleep 5;

rm temp.txt #huge file, get rid of it.


/mnt/sequencedb/PopGen/barbara/simulations/scripts/sample_chr.pl 1 outp1.txt >samplep1;



/mnt/sequencedb/PopGen/barbara/simulations/scripts/sample_chr.pl 60 outp2.txt > samplep2;

/mnt/sequencedb/PopGen/barbara/simulations/scripts/sample_chr.pl 60 outp3.txt > samplep3;

rm outp*

cat samplep1 samplep2 samplep3 > slim_sample.txt

rm samplep*


