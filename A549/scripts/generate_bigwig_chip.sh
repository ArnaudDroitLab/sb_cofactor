MUGQIC_DIR=$1

for i in $MUGQIC_DIR/alignment/*/*.sorted.dup.bam
do
    samplename=`basename $i .sorted.dup.bam`
    if [ ! -e $MUGQIC_DIR/tracks/$samplename.bw ]
    then     
        mkdir -p $MUGQIC_DIR/jobs
        script=$MUGQIC_DIR/jobs/$samplename.make_bigwig.sh
        cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A $RAP_ID
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
module load mugqic/python/2.7.12
bamCoverage -e 200 --binSize 5 -p 16 --normalizeUsingRPKM \
    -b $i \
    -o $MUGQIC_DIR/tracks/$samplename.bw
EOF
        workdir=`pwd`
        qsub $script -o $script.stdout -e $script.stderr -d $workdir
    fi
done
