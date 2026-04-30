# taps_somatic
Somatic SNV Caller for TAPS+ datasets

gcc taps_somatic.c -O3 -lhts -lm -o taps

cfDNA
./taps \
  -t tumor.bam \
  -n normal.bam \
  -r ref.fa \
  --min-depth 200 \
  --min-vaf 0.005 \
  --llr 4.5 \
  --dispersion 30 \
  --contamination 0.02

Tumor
./taps \
  -t tumor.bam \
  -n normal.bam \
  -r ref.fa \
  --min-depth 20 \
  --min-vaf 0.02 \
  --llr 8 \
  --dispersion 80 \
  --contamination 0.005
