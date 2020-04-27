#!/bin/bash

#=================================================================================================================
# ACValidator launcher script
# This script launches ACValidator for all coordinates in ${coordinateFile} for sample ${sample}
#
#
# Example Invocation:
# ./ACV_launcher_v2.sh inputSam coordinateList.txt window_size logFileName.txt
#
#
# Author:
# Shobana Sekar, Liang lab
# Translational Genomics Research Institute, Phoenix, AZ
#
# Date:
# December 2018
#=================================================================================================================

echo -n "Would you like to compute insert size from input sam file? (y/n) "
read answer

if [ "$answer" == "y" ] ;then

    if [ "$#" -ne 3 ]; then
        echo "Incorrect number of arguments, please provide inputSam, coordinateList and logFilename"
        exit 1
    fi

    sample="$1" # Sample to validate the junctions on (provide part of the sample name without the ".sam")
    coordinateFile="$2" # circRNA coordinates to validate in the sample list
    logFile="$3"

    echo "Calculating insert size from input SAM file..."
    PICARDPATH=/packages/picard-tools/1.128/picard.jar
    echo "java -jar ${PICARDPATH} CollectInsertSizeMetrics I=${sample}.sam O=${sample}_insert_size.txt H=${sample}_insert_size.hist.pdf M=0.05"
    java -jar ${PICARDPATH} CollectInsertSizeMetrics I=${sample}.sam O=${sample}_insert_size.txt H=${sample}_insert_size.hist.pdf M=0.05

    insertSize=`grep -A 1 "MEDIAN_INSERT_SIZE" ${sample}_insert_size.txt | tail -n 1 | awk '{print $1}'`

    echo "Suggested window size is: ${insertSize}"

    for line in `cat ${coordinateFile}`
    do
        circRNA=${line}
        echo Processing ${circRNA}
        echo "Calling python ACValidator_v1.py -i ${sample} -c ${circRNA} -w ${insertSize} --log-filename ${logFile}"
        ACValidator -i ${sample} -c ${circRNA} -w ${insertSize} --log-filename ${logFile}
    done

else

    echo No
    if [ "$#" -ne 4 ]; then
        echo "Incorrect number of arguments, please provide inputSam, coordinateList,windowSize and logFilename"
        exit 1
    fi

    sample="$1" # Sample to validate the junctions on (provide part of the sample name without the ".sam")
    coordinateFile="$2" # circRNA coordinates to validate in the sample list
    windowSize="$3" # user-defined window size
    logFile="$4"

    for line in `cat ${coordinateFile}`
    do
	circRNA=${line}
	echo Processing ${circRNA}
	echo "Calling python ACValidator_v1.py -i ${sample} -c ${circRNA} -w ${windowSize} --log-filename ${logFile}"
	ACValidator -i ${sample} -c ${circRNA} -w ${windowSize} --log-filename ${logFile}
    done

fi
