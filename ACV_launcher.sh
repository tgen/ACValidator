#!/bin/bash

#=================================================================================================================
# ACValidator launcher script
# This script launches ACValidator for all coordinates in ${coordinateFile} for sample ${sample}
#
#
# Example Invocation:
# ./ACValidator_launcher.sh inputSam coordinateList.txt window_size logFileName.txt
#
#
# Author:
# Shobana Sekar, Liang lab
# Translational Genomics Research Institute, Phoenix, AZ
#
# Date:
# December 2018
#=================================================================================================================


if [ "$#" -ne 4 ]; then
    echo "Incorrect number of arguments, please provide inputSam, coordinateList and logFilename"
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
	python ACValidator_v1.py -i ${sample} -c ${circRNA} -w ${windowSize} --log-filename ${logFile}
done
