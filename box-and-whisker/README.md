# box-plotv2.0.py
A python2 program for generating box and whisker plots of CDR3 length distributions from a list of CDR3 lengths or what we call a "clonotype file."

## Input type 1 ("clonotype file"):
Examples files are given in `a.dat.gz`, `b.dat.gz`, etc. Each line in the file contains the V and J germline gene family, CDR3, and number of associated somatic variants.

TRBV7-8 TRBJ2-7 ASSFLRYSPYEQY 1695020

## Input type 2 ("data file"):
CDR3 length (one entry per clonotype)
11
13
11
14
8

## Command line options:

### {--clonotype-file | --data-file}

### --files [list of files (e.g., a.dat.gz b.dat.gz)]

### --hex-colors [list of strings in format '#009999' '#ff9999' (arg count must correspond to file count)]

### --y-limits 0 30

### --y-tick-markers 0 5 10 15 20

### --png-file=BoxPlot1.png

### --pdf-file=BoxPlot1.pdf

### --summary-statistics-outfile=statistics.dat.gz \

### --scatter-point-marker=. (optional)

### --scatter-point-size=20 (optional)
