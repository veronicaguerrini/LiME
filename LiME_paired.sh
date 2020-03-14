#!/bin/bash
#usage: LightMetaEbwt_paired.sh File1_F File1_RC File2_F File2_RC output numReads numGenomes readLen threads

#Choose the steps you want to execute
step1=1   #=1 to execute Step 1.
step2=1   #=1 to execute Step 2.
step3=1   #=1 to execute Step 3.

#Default parameters
alpha=16
beta=0.25
TaxRank=1

#Input
FastaFile1_F=$1
FastaFile1_RC=$2
FastaFile2_F=$3
FastaFile2_RC=$4
output=$5
numReads=$6
numGenomes=$7
LineageFile=$8
readLen=$9
threads=${10}


echo "input File1_F: "$FastaFile1_F
echo "input File1_RC: "$FastaFile1_RC
echo "input File2_F: "$FastaFile2_F
echo "input File2_RC: "$FastaFile2_RC
echo "output: "$output
echo "numReads: "$numReads
echo "numGenomes: "$numGenomes
echo "LineageFile: "$LineageFile
echo "readLen: "$readLen
echo "threads: "$threads

deleteFile=1

#First step
if [ $step1 -eq 1 ]
then

nohup /usr/bin/time -v ./ClusterLCP $FastaFile1_F $numReads $numGenomes $alpha $threads > "ClusterLCP_File1_F.stdout" 2> "ClusterLCP_File1_F.stderr" &

nohup /usr/bin/time -v ./ClusterLCP $FastaFile1_RC $numReads $numGenomes $alpha $threads > "ClusterLCP_File1_RC.stdout" 2> "ClusterLCP_File1_RC.stderr" &

nohup /usr/bin/time -v ./ClusterLCP $FastaFile2_F $numReads $numGenomes $alpha $threads > "ClusterLCP_File2_F.stdout" 2> "ClusterLCP_File2_F.stderr" &

nohup /usr/bin/time -v ./ClusterLCP $FastaFile2_RC $numReads $numGenomes $alpha $threads > "ClusterLCP_File2_RC.stdout" 2> "ClusterLCP_File2_RC.stderr" &

#wait the previous processes
wait

fi

#Second Step

if [ $step2 -eq 1 ]
then

nohup /usr/bin/time -v ./ClusterBWT_DA $FastaFile1_F $readLen $beta $threads > "ClusterBWT_DA_File1_F.stdout" 2> "ClusterBWT_DA_File1_F.stderr"

nohup /usr/bin/time -v ./ClusterBWT_DA $FastaFile1_RC $readLen $beta $threads > "ClusterBWT_DA_File1_RC.stdout" 2> "ClusterBWT_DA_File1_RC.stderr"

nohup /usr/bin/time -v ./ClusterBWT_DA $FastaFile2_F $readLen $beta $threads > "ClusterBWT_DA_File2_F.stdout" 2> "ClusterBWT_DA_File2_F.stderr"

nohup /usr/bin/time -v ./ClusterBWT_DA $FastaFile2_RC $readLen $beta $threads > "ClusterBWT_DA_File2_RC.stdout" 2> "ClusterBWT_DA_File2_RC.stderr"

#wait the previous processes
wait

fi

#Third step

if [ $step3 -eq 1 ]
then
/usr/bin/time -v ./Classify 4 $FastaFile1_F".res" $FastaFile1_RC".res" $FastaFile2_F".res" $FastaFile2_RC".res" $numReads $numGenomes $output $LineageFile $TaxRank $threads > "Classify_"$output".stdout" 2> "Classify_"$output".stderr"

fi

#deleteFile
if [ $deleteFile -eq 1 ]
then
rm -f ./Datasets/$FastaFile1_F".res.bin"
rm -f ./Datasets/$FastaFile1_F".res.pos"
rm -f ./Datasets/$FastaFile1_RC".res.bin"
rm -f ./Datasets/$FastaFile1_RC".res.pos"
rm -f ./Datasets/$FastaFile2_F".res.bin"
rm -f ./Datasets/$FastaFile2_F".res.pos"
rm -f ./Datasets/$FastaFile2_RC".res.bin"
rm -f ./Datasets/$FastaFile2_RC".res.pos"

fi

