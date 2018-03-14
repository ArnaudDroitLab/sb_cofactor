OLD_OUTPUT=/home/efournie/projects/def-stbil30/Working_Directory/Eric/A549_Michele/output/chip-pipeline-GRCh38
NEW_OUTPUT=/home/efournie/projects/def-stbil30/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38

OLD_PREFIX=BRD4_NA_EtOH_Rep2
NEW_PREFIX=A549_CTRL_BRD4_rep2
mkdir -p $NEW_OUTPUT/alignment/$NEW_PREFIX
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bai       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bai
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam.flagstat  $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam.flagstat
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.metrics   $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.metrics

OLD_PREFIX=BRD4_NA_Dex_Rep2
NEW_PREFIX=A549_DEX_BRD4_rep2
mkdir -p $NEW_OUTPUT/alignment/$NEW_PREFIX
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bai       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bai
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam.flagstat  $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam.flagstat
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.metrics   $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.metrics

OLD_PREFIX=CDK9_NA_EtOH_Rep2
NEW_PREFIX=A549_CTRL_CDK9_rep2
mkdir -p $NEW_OUTPUT/alignment/$NEW_PREFIX
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bai       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bai
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam.flagstat  $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam.flagstat
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.metrics   $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.metrics

OLD_PREFIX=CDK9_NA_Dex_Rep2
NEW_PREFIX=A549_DEX_CDK9_rep2
mkdir -p $NEW_OUTPUT/alignment/$NEW_PREFIX
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bai       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bai
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam.flagstat  $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam.flagstat
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.metrics   $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.metrics

OLD_PREFIX=WCE_NA_Dex_Rep2
NEW_PREFIX=A549_DEX_WCE_rep2
mkdir -p $NEW_OUTPUT/alignment/$NEW_PREFIX
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bai       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bai
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam.flagstat  $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam.flagstat
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.metrics   $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.metrics

OLD_PREFIX=WCE_NA_EtOH_Rep2
NEW_PREFIX=A549_CTRL_WCE_rep3
mkdir -p $NEW_OUTPUT/alignment/$NEW_PREFIX
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bai       $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bai
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.bam.flagstat  $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.bam.flagstat
mv $OLD_OUTPUT/alignment/$OLD_PREFIX/$OLD_PREFIX.sorted.dup.metrics   $NEW_OUTPUT/alignment/$NEW_PREFIX/$NEW_PREFIX.sorted.dup.metrics



