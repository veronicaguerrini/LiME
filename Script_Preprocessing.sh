#!/bin/bash

########################################################
########################################################
#Script for building the three data structures (ebwt, da and lcp)
#for the entire sequence set of reads and genomes starting
#from a read fasta file and a genome fasta file
########################################################
########################################################
build_DS_reads=1      #=1 to build data structures for the read sets
build_DS_genomes=1    #=1 to build data structures for the genome set
build_merge=1         #=1 to merge data structures
########################################################
########################################################

#Please, modify the following commands according to the fasta file paths and names, and the read collection type (paired-end or single-end)

paired=1 #1=paired-end, 0=single-end
PathReference="./example"
FastaReference="refs.fasta"
PathDataset="./example"
FastaDatasetR1="reads_1.fasta"
FastaDatasetR2="reads_2.fasta"

#Please, set the following variables in case you want to truncate the bwt merge after k iterations
#Recall that using the option --trlcp k, eGap computes an LCP array in which all values greater than k are replaced by the value k.

truncate=0 #1=merge WITH option --trlcp k, 0=merge WITHOUT option --trlcp k
k=20 #set the parameter k

############################
############################
############################
echo "FastaReference: "$PathReference/$FastaReference
#compute the number of references
nRefs=$(grep ">" $PathReference/$FastaReference | wc -l)
echo "Number of genomes: "$nRefs

echo "FastaDatasetR1: "$PathDataset/$FastaDatasetR1
#compute the number of reads
nReads=$(grep ">" $PathDataset/$FastaDatasetR1 | wc -l)
echo "Number of reads: "$nReads

pathseqtk="./Preprocessing/seqtk"
pathEGSA="./Preprocessing/egsa"
pathBCR="./Preprocessing/BCR_LCP_GSA"
pathEGAP="./Preprocessing/egap"

#reverse complement R1
FastaDatasetR1RC="$(basename "$FastaDatasetR1" .fasta)_RC.fasta"
echo "FastaDatasetR1_RC: "$PathDataset/$FastaDatasetR1RC

$pathseqtk/seqtk seq -r "$PathDataset/$FastaDatasetR1" > "$PathDataset/$FastaDatasetR1RC"
#####

#reverse complement R2
if [ $paired -eq 1 ]
then
    echo -e "\nPaired-end collection\nFastaDatasetR2: "$PathDataset/$FastaDatasetR2
    FastaDatasetR2RC="$(basename "$FastaDatasetR2" .fasta)_RC.fasta"
    echo "FastaDatasetR2_RC: "$PathDataset/$FastaDatasetR2RC

    $pathseqtk/seqtk seq -r "$PathDataset/$FastaDatasetR2" > "$PathDataset/$FastaDatasetR2RC"
fi
#####

g++ create_docs.cpp -o ./Preprocessing/create_docs

#BCR for eBWT/DA for reads (short sequences)

if [ $build_DS_reads -eq 1 ]
then
    mkdir DS_reads

    echo -e "\nComputing eBWT/DA for set R1"
    echo "pathBCR: "$pathBCR
    echo "Start BCR..."

    /usr/bin/time -v $pathBCR/BCR_LCP_GSA $PathDataset/$FastaDatasetR1 ./DS_reads/$FastaDatasetR1 > "./DS_reads/BCR_$(basename "$FastaDatasetR1" .fasta).stdout" 2> "./DS_reads/BCR_$(basename "$FastaDatasetR1" .fasta).stderr"

    ./Preprocessing/create_docs $PathDataset/$FastaDatasetR1 $nReads
    mv $PathDataset/$FastaDatasetR1".docs" ./DS_reads
    mv ./DS_reads/$FastaDatasetR1".da" ./DS_reads/$FastaDatasetR1".4.da"

    echo -e "\nComputing eBWT/DA for set R1_RC"
    echo "Start BCR..."

    /usr/bin/time -v $pathBCR/BCR_LCP_GSA $PathDataset/$FastaDatasetR1RC ./DS_reads/$FastaDatasetR1RC > "./DS_reads/BCR_$(basename "$FastaDatasetR1RC" .fasta).stdout" 2> "./DS_reads/BCR_$(basename "$FastaDatasetR1RC" .fasta).stderr"

    ./Preprocessing/create_docs $PathDataset/$FastaDatasetR1RC $nReads
    mv $PathDataset/$FastaDatasetR1RC".docs" ./DS_reads
    mv ./DS_reads/$FastaDatasetR1RC".da" ./DS_reads/$FastaDatasetR1RC".4.da"

    if [ $paired -eq 1 ]
    then
        echo -e "\nComputing eBWT/DA for set R2"
        echo "Start BCR..."

        /usr/bin/time -v $pathBCR/BCR_LCP_GSA $PathDataset/$FastaDatasetR2 ./DS_reads/$FastaDatasetR2 > "./DS_reads/BCR_$(basename "$FastaDatasetR2" .fasta).stdout" 2> "./DS_reads/BCR_$(basename "$FastaDatasetR2" .fasta).stderr"

        ./Preprocessing/create_docs $PathDataset/$FastaDatasetR2 $nReads
        mv $PathDataset/$FastaDatasetR2".docs" ./DS_reads
        mv ./DS_reads/$FastaDatasetR2".da" ./DS_reads/$FastaDatasetR2".4.da"

        echo -e "\nComputing eBWT/DA for set R2_RC"
        echo "Start BCR..."

        /usr/bin/time -v $pathBCR/BCR_LCP_GSA $PathDataset/$FastaDatasetR2RC ./DS_reads/$FastaDatasetR2RC > "./DS_reads/BCR_$(basename "$FastaDatasetR2RC" .fasta).stdout" 2> "./DS_reads/BCR_$(basename "$FastaDatasetR2RC" .fasta).stderr"

        ./Preprocessing/create_docs $PathDataset/$FastaDatasetR2RC $nReads
        mv $PathDataset/$FastaDatasetR2RC".docs" ./DS_reads
        mv ./DS_reads/$FastaDatasetR2RC".da" ./DS_reads/$FastaDatasetR2RC".4.da"
    fi
    rm ./DS_reads/*.len
    rm ./DS_reads/*.info
fi
####
####

#EGSA for eBWT/DA for genomes (long sequences)

if [ $build_DS_reads -eq 1 ]
then
    mkdir DS_genomes

    echo -e "\nComputing eBWT/DA for genomes"
    echo "pathEGSA: "$pathEGSA
    echo "Start EGSA..."

    /usr/bin/time -v $pathEGSA/egsa $PathReference/$FastaReference $nRefs > "./DS_genomes/EGSA_$(basename "$FastaReference" .fasta).stdout" 2> "./DS_genomes/EGSA_$(basename "$FastaReference" .fasta).stderr"

    rm -fr $PathReference/partition
    rm -fr $PathReference/tmp

    ./Preprocessing/create_docs $PathReference/$FastaReference $nRefs
    echo "Start EGSAtoBCR..."
    ./EGSAtoBCR $PathReference/$FastaReference $nRefs
    mv $PathReference/$FastaReference".ebwt" ./DS_genomes
    mv $PathReference/$FastaReference".da" ./DS_genomes/$FastaReference".4.da"
    mv $PathReference/$FastaReference".docs" ./DS_genomes
    rm $PathReference/$FastaReference".lcp"
    rm $PathReference/$FastaReference"."$nRefs".gesa"
fi
####
####

#EGAP to merge BWT files and compute LCP and DA of the entire collection
if [ $build_DS_merge -eq 1 ]
then
    mkdir DS_merge
    echo -e "\nComputing eBWT/DA/LCP for set reads+genomes"
    echo "pathEGAP: "$pathEGAP
    echo "Start EGAP..."
    if [ $truncate -eq 1 ]
    then
        outmergeF="$(basename "$FastaDatasetR1" .fasta)+Refs_tr"$k".fasta"
        /usr/bin/time -v $pathEGAP/eGap -m 4096 --bwt -o ./DS_merge/$outmergeF ./DS_reads/$FastaDatasetR1".ebwt" ./DS_genome/$FastaReference".ebwt" --lbytes 4 --da --dbytes 4 --docs --trlcp $k > "eGap_merge_$(basename "$outmergeF" .fasta).stdout" 2> "eGap_merge_$(basename "$outmergeF" .fasta).stderr"
        outmergeRC="$(basename "$FastaDatasetR1RC" .fasta)+Refs_tr"$k".fasta"
        /usr/bin/time -v $pathEGAP/eGap -m 4096 --bwt -o ./DS_merge/$outmergeRC ./DS_reads/$FastaDatasetR1RC".ebwt" ./DS_genome/$FastaReference".ebwt" --lbytes 4 --da --dbytes 4 --docs --trlcp $k > "eGap_merge_$(basename "$outmergeRC" .fasta).stdout" 2> "eGap_merge_$(basename "$outmergeRC" .fasta).stderr"
    else
        outmergeF="$(basename "$FastaDatasetR1" .fasta)+Refs.fasta"
        /usr/bin/time -v $pathEGAP/eGap -m 4096 --bwt -o ./DS_merge/$outmergeF ./DS_reads/$FastaDatasetR1".ebwt" ./DS_genome/$FastaReference".ebwt" --lcp --lbytes 4 --da --dbytes 4 --docs > "eGap_merge_$(basename "$outmergeF" .fasta).stdout" 2> "eGap_merge_$(basename "$outmergeF" .fasta).stderr"
        outmergeRC="$(basename "$FastaDatasetR1RC" .fasta)+Refs.fasta"
        /usr/bin/time -v $pathEGAP/eGap -m 4096 --bwt -o ./DS_merge/$outmergeRC ./DS_reads/$FastaDatasetR1RC".ebwt" ./DS_genome/$FastaReference".ebwt" --lcp --lbytes 4 --da --dbytes 4 --docs > "eGap_merge_$(basename "$outmergeRC" .fasta).stdout" 2> "eGap_merge_$(basename "$outmergeRC" .fasta).stderr"
    fi
    mv ./DS_merge/$outmergeF".4.lcp" ./DS_merge/$outmergeF".lcp"
    mv ./DS_merge/$outmergeRC".4.lcp" ./DS_merge/$outmergeRC".lcp"
    mv ./DS_merge/$outmergeF".4.da" ./DS_merge/$outmergeF".da"
    mv ./DS_merge/$outmergeRC".4.da" ./DS_merge/$outmergeRC".da"
    mv ./DS_merge/$outmergeF".bwt" ./DS_merge/$outmergeF".ebwt"
    mv ./DS_merge/$outmergeRC".bwt" ./DS_merge/$outmergeRC".ebwt"

    if [ $paired -eq 1 ]
    then
        if [ $truncate -eq 1 ]
        then
            outmergeF="$(basename "$FastaDatasetR2" .fasta)+Refs_tr"$k".fasta"
            /usr/bin/time -v $pathEGAP/eGap -m 4096 --bwt -o ./DS_merge/$outmergeF ./DS_reads/$FastaDatasetR2".ebwt" ./DS_genome/$FastaReference".ebwt" --lbytes 4 --da --dbytes 4 --docs --trlcp $k > "eGap_merge_$(basename "$outmergeF" .fasta).stdout" 2> "eGap_merge_$(basename "$outmergeF" .fasta).stderr"
            outmergeRC="$(basename "$FastaDatasetR2RC" .fasta)+Refs_tr"$k".fasta"
            /usr/bin/time -v $pathEGAP/eGap -m 4096 --bwt -o ./DS_merge/$outmergeRC ./DS_reads/$FastaDatasetR2RC".ebwt" ./DS_genome/$FastaReference".ebwt" --lbytes 4 --da --dbytes 4 --docs --trlcp $k > "eGap_merge_$(basename "$outmergeRC" .fasta).stdout" 2> "eGap_merge_$(basename "$outmergeRC" .fasta).stderr"

        else
            outmergeF="$(basename "$FastaDatasetR2" .fasta)+Refs.fasta"
            /usr/bin/time -v $pathEGAP/eGap -m 4096 --bwt -o ./DS_merge/$outmergeF ./DS_reads/$FastaDatasetR2".ebwt" ./DS_genome/$FastaReference".ebwt" --lcp --lbytes 4 --da --dbytes 4 --docs > "eGap_merge_$(basename "$outmergeF" .fasta).stdout" 2> "eGap_merge_$(basename "$outmergeF" .fasta).stderr"
            outmergeRC="$(basename "$FastaDatasetR2RC" .fasta)+Refs.fasta"
            /usr/bin/time -v $pathEGAP/eGap -m 4096 --bwt -o ./DS_merge/$outmergeRC ./DS_reads/$FastaDatasetR2RC".ebwt" ./DS_genome/$FastaReference".ebwt" --lcp --lbytes 4 --da --dbytes 4 --docs > "eGap_merge_$(basename "$outmergeRC" .fasta).stdout" 2> "eGap_merge_$(basename "$outmergeRC" .fasta).stderr"
        fi
        mv ./DS_merge/$outmergeF".4.lcp" ./DS_merge/$outmergeF".lcp"
        mv ./DS_merge/$outmergeRC".4.lcp" ./DS_merge/$outmergeRC".lcp"
        mv ./DS_merge/$outmergeF".4.da" ./DS_merge/$outmergeF".da"
        mv ./DS_merge/$outmergeRC".4.da" ./DS_merge/$outmergeRC".da"
        mv ./DS_merge/$outmergeF".bwt" ./DS_merge/$outmergeF".ebwt"
        mv ./DS_merge/$outmergeRC".bwt" ./DS_merge/$outmergeRC".ebwt"
    fi
    rm ./DS_merge/*.docs
fi
echo "Done."
####
####

