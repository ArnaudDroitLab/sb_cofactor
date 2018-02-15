for nr in LCC9
do
    for region in Promoter First Other Intergenic
    do
        Rscript scripts/CIHR-metagene-command.R $nr $region &
    done
done

