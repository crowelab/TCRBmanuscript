# check-for-overlapsv0.1.py
A python2 script for determining shared clonotypes between datasets

## Clonotype file definition:
Examples files are given in TCR-HIP directory. Each line in the file contains the V and J germline gene family, CDR3, and number of associated somatic variants.

TRBV7-8 TRBJ2-7 ASSFLRYSPYEQY 1695020

## Command line options:
### --dump-clonotypes
Creates a list of clonotypes shared between the input sets

### --outfile=statistics.dat.gz
Save summary of overlapping clonotypes

### --common-clonotypes-file=commons.dat.gz
Write shared clonotypes to file

### --clonotype-files A.dat.gz B.dat.gz C.dat.gz
Input data files
