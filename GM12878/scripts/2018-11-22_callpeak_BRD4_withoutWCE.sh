# GEO GSE62912
# Callpeak on SRR1636861_BRD4
# without WCE (not available)

mkdir -p output/chip-pipeline-GRCh38/peak_call/SRR1636861_BRD4

macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  output/chip-pipeline-GRCh38/alignment/SRR1636861_BRD4/SRR1636861_BRD4.sorted.dup.bam \
  --name output/chip-pipeline-GRCh38/peak_call/SRR1636861_BRD4/SRR1636861_BRD4 \
  >& output/chip-pipeline-GRCh38/peak_call/SRR1636861_BRD4/SRR1636861_BRD4.diag.macs.out
