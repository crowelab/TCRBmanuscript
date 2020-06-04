
#export OMP_NUM_THREADS=1
MAIN=`pwd`
BINARY=paired-subsampling.py

python2 $BINARY --clonotype-file1 $MAIN/TestData/AB-HELIX-HIP1-TCRB.dat.gz \
               --clonotype-file2 $MAIN/TestData/ADAPTIVE-HIP1-TCRB.dat.gz \
               --beta --beta-metrics=percentage-overlap \
	       --standardize-sample-size --initialize=10 --number-sweeps=1000 \
               --outfile=raw-statistics-HIP1.dat.gz --outfile-summary-statistics=summary-statistics-HIP1.dat.gz 

