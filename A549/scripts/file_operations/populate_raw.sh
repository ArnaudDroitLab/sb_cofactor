for f in /gs/project/fhq-091-aa/Raw_Data/RNA_2014-08-19/*.bam
do
    newname=`basename $f | sed -e 's/HI.2031.00..Index_.*SB_A549-//'`
    ln -s $f raw/$newname
done
