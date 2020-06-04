HIP=/SourceCode/Python/Plots/COLLECTIONS/HIPS
BINARY=morisita-index-vjV3.0.py
python2 $BINARY --clonotype-files1 $HIP/HIP1/AB-HELIX-HIP1-TCRB.dat.gz \
               --clonotype-files2 $HIP/HIP1/ADAPTIVE-HIP1-TCRB.dat.gz \
               --vj-frequency-only --create-plot \
               --hex-colors '#009999' 
#SHOULD BE 0.25
python2 $BINARY --clonotype-files1 $HIP/HIP2/AB-HELIX-HIP2-TCRB.dat.gz \
               --clonotype-files2 $HIP/HIP2/ADAPTIVE-HIP2-TCRB.dat.gz \
               --vj-frequency-only --create-plot \
               --hex-colors '#ff9933' 
#should be 0.71

python2 $BINARY --clonotype-files1 $HIP/HIP3/AB-HELIX-HIP3-TCRB.dat.gz \
               --clonotype-files2 $HIP/HIP3/ADAPTIVE-HIP3-TCRB.dat.gz \
               --vj-frequency-only --create-plot \
               --hex-colors '#ff44ab'
#should be 0.62

