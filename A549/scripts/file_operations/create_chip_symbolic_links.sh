#!/bin/bash

# This script creates the symbolic link from the sequencing results
# in def-stbil30/Raw_Data to the raw directory within the cofactor git
# tree.

# The expected directory structure is
# STEVE_ROOT
#  |-> Raw_Data
#  |-> Working_Directory
#       |-> Eric
#            |-> CofactorGIT
#             


pushd raw/chip-seq/
ln -s ../../../../../../Raw_Data/ChIP_2014-02-04/*A549*.bam ./
ln -s ../../../../../../Raw_Data/ChIP_2014-06-10/*A549*.bam ./
ln -s ../../../../../../Raw_Data/ChIP_2014-12-22/*A549*.bam ./
ln -s ../../../../../../Raw_Data/ChIP-2017-12-06/*A549*CDK9*.fastq.gz ./
ln -s ../../../../../../Raw_Data/ChIP-2017-12-06/*A549*BRD4*.fastq.gz ./
ln -s ../../../../../../Raw_Data/ChIP-2017-12-06/*A549*WCE_NA*.fastq.gz ./
popd