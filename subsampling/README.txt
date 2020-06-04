#Subsampling script is used to standardize sample sizes before checking
#for percentage of overlapping V3J clonotypes. The files provided in this
#directory are for testing purposes. Please make changes to the Bash 
#script below for use:
# --clonotype-file1                                           : list of clonotype files to compare (must be ordered with clonotype-file2 flag)
# --clonotype-file1                                           : list of clonotype files to compare (must be ordered with clonotype-file1 flag)
# --standardize-sample-size                                   : required flag in order to standardize the sizes of both files
# --beta-metrics=percentage-overlap                           : required to perform percentage overlap after each random shuffles
# --intialize=10                                              : will perform 10 random shuffles before perform full set of shuffles
# --number-sweeps=1000                                        : will perform 1000 shuffles each time calculating the percentage overlaps
# --outfile=raw-statistics-HIP1.dat.gz                        : dump out percentage overlaps after each random shuffle into file
# --outfile-summary-statistics=summary-statistics-HIP1.dat.gz : dump out mean, median and standard deviation of all shuffling to file

MAIN=`pwd`
BINARY=paired-subsampling.py

python2 $BINARY --clonotype-file1 $MAIN/TestData/AB-HELIX-HIP1-TCRB.dat.gz \
               --clonotype-file2 $MAIN/TestData/ADAPTIVE-HIP1-TCRB.dat.gz \
               --beta --beta-metrics=percentage-overlap \
	       --standardize-sample-size --initialize=10 --number-sweeps=1000 \
               --outfile=raw-statistics-HIP1.dat.gz --outfile-summary-statistics=summary-statistics-HIP1.dat.gz 

