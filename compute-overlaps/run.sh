MAIN=`pwd`
COLLECTIONS=$MAIN/TCR-HIP/
BINARY=check-for-overlapsv0.1.py
python2 $BINARY --dump-clonotypes \
               --outfile=statistics.dat.gz \
               --common-clonotypes-file=commons.dat.gz \
               --clonotype-files $COLLECTIONS/A.dat.gz \
	                         $COLLECTIONS/B.dat.gz \
                                 $COLLECTIONS/C.dat.gz 
 
