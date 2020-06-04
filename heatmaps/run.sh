BASE=/home/sotocs/SourceCode/Python/Plots/COLLECTIONS/SUBSETS

#################################
# SCRIPT FOR CREATING HEATMAPS  #
#################################
python2 heatmap-generation.py  --clonotype-files $BASE/CD8+NAIVE/HIP2-NAIVE-CD8+.dat.gz \
                                                $BASE/CD8+NAIVE/HIP3-NAIVE-CD8+.dat.gz \
                                                --dump-csv \
                                                --output-csv-files HIP2-NAIVE-CD8+.csv.gz \
                                                                   HIP3-NAIVE-CD8+.csv.gz \
                                                --output-aggregate-csv-file=aggregate.csv.gz \
                                                --heatmap-output-pdf=test.pdf \
                                                --heatmap-output-png=test.png
