cd raw/ENCODE-chip

for i in `sed 1d readset.txt | cut -f 7`
do
    wget $i
done
