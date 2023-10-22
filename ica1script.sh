#!/usr/bin/bash

# 0 copy the ICA1 files in my working directory
cp -r /localdisk/data/BPSM/ICA1 ${PWD}
cd ${PWD}/ICA1

# 1 a quality check using fastqc

# clear previous output files
rm -fr ${PWD}/qcoutput
# creat an output folder for fastqc
mkdir ${PWD}/qcoutput

# the input data is FASTQ format, so use "-f fastq" ; not extracts the data files; outputfile is 'outputdir'
fastqc -f fastq -noextract -t 8 -o ${PWD}/qcoutput ${PWD}/fastq/*.gz

# 2 assess

# 3 bowtie2 mapping; convert outputfile to bam format and sort

# 3.1 uncompress 'TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz' file
rm -f ${PWD}/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta
gunzip -k ${PWD}/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz

# 3.2 build Bowtie2 index
rm -f ${PWD}/Tcongo_genome/*.bt2
bowtie2-build ${PWD}/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta ${PWD}/Tcongo_genome/TcongolenseIL3000_2019_Genome

# 3.3 use bowtie2 to align each pair of sequences; convert outputfile to BAM format and sort BAM files

rm -f ${PWD}/fastq/sequence_file.txt
# extract the file name of the sequence that requires alignment and save it to sequence_file.txt
ls ${PWD}/fastq | grep '.gz' > ${PWD}/fastq/sequence_file.txt

IFS=$'_';
rm -fr ${PWD}/alignoutput
mkdir ${PWD}/alignoutput
alignoutput=${PWD}/alignoutput
while read Tcopart1 Tcopart2
do
bowtie2 -p 10 -x ${PWD}/Tcongo_genome/TcongolenseIL3000_2019_Genome -1 ${PWD}/fastq/${Tcopart1}_1.fq.gz  -2 ${PWD}/fastq/${Tcopart1}_2.fq.gz -S ${alignoutput}/${Tcopart1}.sam
samtools view -bh -S ${alignoutput}/${Tcopart1}.sam > ${alignoutput}/${Tcopart1}.bam
samtools sort -O BAM -m 8g -@ 6 -o ${alignoutput}/${Tcopart1}_sort.bam ${alignoutput}/${Tcopart1}.bam
done < ${PWD}/fastq/sequence_file.txt
unset IFS

# 4 generate counts data using bedtools
IFS=$'_'
rm -fr ${PWD}/genecounts
mkdir ${PWD}/genecounts
genecounts=${PWD}/genecounts
while read Tcopart1 Tcopart2
do
bedtools coverage -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b ${PWD}/fastq/alignoutput/${Tcopart1}_sort.bam > ${genecounts}/${Tcopart1}.gene.counts.txt
done < ${PWD}/fastq/sequence_file.txt
unset IFS

# 5 calculate average of the counts

# 5.1 group the samples

rm -fr ${PWD}/group

# create a folder for grouping
mkdir ${PWD}/group
group=${PWD}/group

# create three arrays for loop generation of grouping files
sampletype=("Clone1" "Clone2" "WT")
sampletime=("0" "24" "48")
treatment=("Uninduced" "Induced")
for item1 in "${sampletype[@]}"
do
    for item2 in "${sampletime[@]}"
    do
		for item3 in "${treatment[@]}"
		do
		cat ${PWD}/fastq/Tco2.fqfiles |grep -w "$item1" | grep -w "$item2" | grep -w "$item3" | cut -f1 >> ${group}/"$item1"_"$item2"_"$item3"
		done
	done
done

# delete empty files
find ${group}/ -type f -empty -delete

# Create a file that records the results of the grouping
for item in $(ls ${group})
do
echo "$item" >> ${group}/groupinfo1.txt
paste -s ${group}/$item >> ${group}/groupinfo2.txt
done

paste ${group}/groupinfo1.txt ${group}/groupinfo2.txt > ${group}/groupinfo3.txt
awk '{gsub("Tco", "Tco-"); print}' ${group}/groupinfo3.txt > ${group}/groupinfo.txt

rm -f ${group}/groupinfo1.txt ${group}/groupinfo2.txt ${group}/groupinfo3.txt

# 5.2 create gene descriptions file
echo -e "name\tgene_description" > ${genecounts}/gene_descriptions.txt
cut -f4,5 ${PWD}/TriTrypDB-46_TcongolenseIL3000_2019.bed >> ${genecounts}/gene_descriptions.txt

# 5.3 calculate the mean expression levels for each group and add gene descriptions

for item in $(ls ${genecounts} | grep "txt")
do
cut -f6 ${genecounts}/$item > ${genecounts}/cut.$item
done

# clear the previous results
rm -f ${genecounts}/*ed*


# the main purpose for this loop: 
# 1) extract the sixth column of each bedtools output file 
# 2) paste them together according to grouping outcomes 
# 3) calculate the mean expression levels
# 4) add gene descriptions

IFS=$'\t';
while read groupname item1 item2 item3 item4
do
# based on the sample items of each group, there are two conditions : 3 items and 4 items
# if there are 3 items
if test -z $item4;
then
# create a header of each group counts file
echo -e "$item1\t$item2\t$item3" > ${genecounts}/${groupname}.temp1
# according to the group, paste the gene counts files together
paste ${genecounts}/cut.$item1.gene.counts.txt ${genecounts}/cut.$item2.gene.counts.txt ${genecounts}/cut.$item3.gene.counts.txt >> ${genecounts}/${groupname}.temp1
# calculate the average of counts for each group
cat ${genecounts}/${groupname}.temp1 | source ${genecounts}/average.sh > ${genecounts}/${groupname}.temp2
# Merge the gene description file and the average of the expression levels file
paste ${genecounts}/gene_descriptions.txt ${genecounts}/${groupname}.temp2 > ${genecounts}/${groupname}.counts

# if there are 4 items
else
echo -e "$item1\t$item2\t$item3\t$item4" > ${genecounts}/${groupname}.temp1
paste ${genecounts}/cut.$item1.gene.counts.txt ${genecounts}/cut.$item2.gene.counts.txt ${genecounts}/cut.$item3.gene.counts.txt ${genecounts}/cut.$item4.gene.counts.txt >> ${genecounts}/${groupname}.temp1
cat ${genecounts}/${groupname}.temp1 | source ${genecounts}/average.sh > ${genecounts}/${groupname}.temp2
paste ${genecounts}/gene_descriptions.txt ${genecounts}/${groupname}.temp2 > ${genecounts}/${groupname}.counts
fi
done < ${group}/groupinfo.txt

unset IFS

rm -f ${genecounts}/cut* ${genecounts}/temp*