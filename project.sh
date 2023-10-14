#!/bin/bash
rm -f ${PWD}/fastq/sequence_file.txt
# extract the file name of the sequence that requires fastqc and save it to sequence_file.txt
ls ${PWD}/fastq | grep '.gz' > ${PWD}/fastq/sequence_file.txt

# a quality check
rm -fr ${PWD}/fastq/outputdir
mkdir ${PWD}/fastq/outputdir
while read Tco
do
# the input data is FASTQ format, so use "-f fastq" ; not extracts the data files; outputfile is 'outputdir'
fastqc -f fastq -noextract -t 8 -o ${PWD}/fastq/outputdir ${PWD}/fastq/${Tco}
done < ${PWD}/fastq/sequence_file.txt

# assess the numbers and quality of the raw sequence data

# uncompress 'TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz' file
rm -f ${PWD}/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta
gunzip -k ${PWD}/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz

# build Bowtie2 index
rm -f ${PWD}/Tcongo_genome/*.bt2
bowtie2-build ${PWD}/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta ${PWD}/Tcongo_genome/TcongolenseIL3000_2019_Genome

