# heatmap-from-clonotypev2.0.py
A python2 script for generating normalized heatmaps in CSV or image (PDF, PNG) format from counts of V and J germline gene usage in an AIRR repertoire.

## Command line options

### --clonotype-files
Example: clonotypes.dat.gz clonotypes2.dat.gz

### --dump-csv
Creates CSV file as well as images

### --output-csv-files 
File names to use with CSV dump
Example: file1.csv.gz file2.csv.gz

### --v-chain-type=IGHV
`['IGHV', 'TRBV']`

### --j-chain-type=IGHJ
`['IGHJ', 'TRBJ']`

### --heatmap-output-pdf=test.pdf
Save PDF file

### --heatmap-output-png=test.png
Save PNG file

### --output-aggregate-csv-file=aggregate.csv.gz
aggregate csv file containing normalized data
