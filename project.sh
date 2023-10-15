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

#bowtie2 mapping and convert outputfile to bam format and sort
IFS="_";
rm -fr ${PWD}/fastq/alignoutput
mkdir ${PWD}/fastq/alignoutput
alignoutput=${PWD}/fastq/alignoutput
while read Tcopart1 Tcopart2
do
bowtie2 -p 10 -x ${PWD}/Tcongo_genome/TcongolenseIL3000_2019_Genome -1 ${PWD}/fastq/${Tcopart1}_1.fq.gz  -2 ${PWD}/fastq/${Tcopart1}_2.fq.gz -S ${alignoutput}/${Tcopart1}.sam
samtools view -bh -S ${alignoutput}/${Tcopart1}.sam > ${alignoutput}/${Tcopart1}.bam
samtools sort -O BAM -m 8g -@ 6 -o ${alignoutput}/${Tcopart1}_sort.bam ${alignoutput}/${Tcopart1}.bam
done < ${PWD}/fastq/sequence_file.txt

rm -rf ${PWD}/genecount
mkdir ${PWD}/genecount

# generate counts data using bedtools
while read Tcopart1 Tcopart2
do
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b ${PWD}/fastq/alignoutput/${Tcopart1}_sort.bam -c > ${PWD}/genecount/${Tcopart1}.gene.counts.txt
done < ${PWD}/fastq/sequence_file.txt
unset IFS


