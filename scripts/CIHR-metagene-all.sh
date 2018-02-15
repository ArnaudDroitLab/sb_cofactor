for nr in ER GR LCC9
do
    for region in Promoter_DE_Up Promoter_DE_Down
    do
        Rscript scripts/CIHR-metagene-command.R $nr $region &
    done
done

