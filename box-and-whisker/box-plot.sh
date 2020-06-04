WORKING=`pwd`
BINARY=$WORKING/box-plotv2.0.py
COLLECTIONS=$WORKING/TCR-HIP/


#######################################
# USED REDUCED FILES FOR TESTING ONLY #
########################################
#zcat $COLLECTIONS/HIPS/HIP1/AB-HELIX-HIP1-TCRB.dat.gz | head -n 1000 > a.dat
#zcat $COLLECTIONS/HIPS/HIP2/AB-HELIX-HIP2-TCRB.dat.gz | head -n 1000 > b.dat
#zcat $COLLECTIONS/HIPS/HIP3/AB-HELIX-HIP3-TCRB.dat.gz | head -n 1000 > c.dat
#zcat $COLLECTIONS/HIPS/HIP1/ADAPTIVE-HIP1-TCRB.dat.gz | head -n 1000 > d.dat
#zcat $COLLECTIONS/HIPS/HIP2/ADAPTIVE-HIP2-TCRB.dat.gz | head -n 1000 > e.dat
#zcat $COLLECTIONS/HIPS/HIP3/ADAPTIVE-HIP3-TCRB.dat.gz | head -n 1000 > f.dat

#pigz -f ?.dat

python2 $BINARY --clonotype-file --files a.dat.gz b.dat.gz c.dat.gz d.dat.gz e.dat.gz f.dat.gz \
--hex-colors '#009999' '#ff9933' '#ff44ab'  '#009999' '#ff9933' '#ff44ab' \
--y-limits 0 30 \
--y-tick-markers 0 5 10 15 20 \
--png-file=BoxPlot1.png --pdf-file=BoxPlot1.pdf


#######################################
# USED REDUCED FILES FOR TESTING ONLY #
########################################
#zcat $COLLECTIONS/HIPS/HIP1/AB-HELIX-HIP1-TCRB.dat.gz | head -n 1000 | gawk '{print length($3)}'  > hip1-ab
#zcat $COLLECTIONS/HIPS/HIP2/AB-HELIX-HIP2-TCRB.dat.gz | head -n 1000 | gawk '{print length($3)}'  > hip2-ab
#zcat $COLLECTIONS/HIPS/HIP3/AB-HELIX-HIP3-TCRB.dat.gz | head -n 1000 | gawk '{print length($3)}'  > hip3-ab
#zcat $COLLECTIONS/HIPS/HIP1/ADAPTIVE-HIP1-TCRB.dat.gz | head -n 1000 | gawk '{print length($3)}'  > hip1-ad
#zcat $COLLECTIONS/HIPS/HIP2/ADAPTIVE-HIP2-TCRB.dat.gz | head -n 1000 | gawk '{print length($3)}'  > hip2-ad
#zcat $COLLECTIONS/HIPS/HIP3/ADAPTIVE-HIP3-TCRB.dat.gz | head -n 1000 | gawk '{print length($3)}'  > hip3-ad


python2 $BINARY  --data-file --files hip1-ab hip2-ab hip3-ab hip1-ad hip2-ad hip3-ad \
	        --hex-colors '#009999' '#ff9933' '#ff44ab' '#009999' '#ff9933' '#ff44ab' \
	        --scatter --y-tick-markers 0 5 10 15 20 25 --y-limits 0 30 --summary-statistics-outfile=statistics.dat.gz \
	        --scatter-point-marker=. --scatter-point-size=20 --png-file=BoxPlot2.png --pdf-file=BoxPlot2.pdf
