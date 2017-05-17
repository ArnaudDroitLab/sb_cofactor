#!/bin/bash

DIR=$1
GENOME=$2

pushd $DIR
cat <<EOF > hub.txt
hub hub_name
shortLabel hub_short_label
longLabel hub_long_label
genomesFile genomes.txt
email Fournier.Eric.2@crchudequebec.ulaval.ca
descriptionUrl 
EOF

echo "genome $GENOME" > genomes.txt
echo "trackDb $GENOME/trackDb.txt" >> genomes.txt

rm $GENOME/trackDb.txt

mkdir $GENOME
cd $GENOME

for track in *.bw
do
    if [ -f $track ]
    then
        bn=`basename $track .bw`
        cat <<EOF >> trackDb.txt
track $bn
bigDataUrl $track
shortLabel $bn
longLabel $bn
type bigWig
autoScale on
alwaysZero on
visibility full

EOF
    fi
done

for track in *.bb
do
    if [ -f $track ]
    then
        bn=`basename $track .bb`
        cat <<EOF >> trackDb.txt
track $bn
bigDataUrl $track
shortLabel $bn
longLabel $bn
type bigBed
visibility dense

EOF
    fi
done

for track in *.bam
do
    if [ -f $track ]
    then
        bn=`basename $track .bam`
        cat <<EOF >> trackDb.txt
track $bn
bigDataUrl $track
shortLabel $bn
longLabel $bn
type bam
visibility dense

EOF
    fi
done


popd

