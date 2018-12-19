# GEO GSE62912
# Callpeak on SRR1636861_BRD4
# without WCE (not available)

mkdir -p output/chip-pipeline-GRCh38/peak_call/SRR1636861_BRD4_withOurWCErep1

macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  output/chip-pipeline-GRCh38/alignment/SRR1636861_BRD4/SRR1636861_BRD4.sorted.dup.bam \
  --control output/chip-pipeline-GRCh38/alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.bam \
  --name output/chip-pipeline-GRCh38/peak_call/SRR1636861_BRD4_withOurWCErep1/SRR1636861_BRD4_withOurWCErep1 \
  >& output/chip-pipeline-GRCh38/peak_call/SRR1636861_BRD4_withOurWCErep1/SRR1636861_BRD4_withOurWCErep1.diag.macs.out
