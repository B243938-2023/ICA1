#!/bin/bash
rm -f sequence_file.txt
# extract the file name of the sequence that requires fastqc and save it to sequence_file.txt
ls ~/ICA1/fastq | grep '.gz' > sequence_file.txt

# a quality check
while read Tco
do
# the input data is FASTQ format,so use '-f fastq';not extracts the data files;outputdir is 'outpufdir'
fastqc -f fastq -noextract -t 8 -o outputdir ${Tco}
done < sequence_file.txt
