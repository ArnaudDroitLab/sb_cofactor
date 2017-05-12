for i in output/chip-pipeline/alignment/*/*.sorted.dup.bam
do
    samplename=`basename $i .sorted.dup.bam`
    mkdir -p output/chip-pipeline/jobs
    script=output/chip-pipeline/jobs/$samplename.make_bigwig.sh
    cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A fhq-091-aa
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
module load mugqic/python/2.7.12
bamCoverage -e 200 --binSize 5 -p 16 --normalizeUsingRPKM \
    -b $i \
    -o output/chip-pipeline/tracks/$samplename.bw
EOF
    workdir=`pwd`
    qsub $script -o $script.stdout -e $script.stderr -d $workdir
done
