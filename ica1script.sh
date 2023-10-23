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

# 2 assess the qc outcomes

#uncompress the fastqc output
unzip -q -d ${PWD}/qcoutput/ './qcoutput/*.zip'

qcoutcomes=$(ls -d ${PWD}/qcoutput/*fastqc)

for qcoutcome in ${qcoutcomes};
	do
 	{
	echo 'quality control results: '${qcoutcome}' ';
	grep -i 'filename' ${qcoutcome}/fastqc_data.txt;
	echo "FAIL :";
	grep 'fail' ${qcoutcome}/fastqc_data.txt;
	echo "WARN :";
	grep 'warn' ${qcoutcome}/fastqc_data.txt;
	} >> ${PWD}/qcoutput/fastqc_summary.txt
done


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

rm -fr ${PWD}/genecounts
mkdir ${PWD}/genecounts
genecounts="${PWD}/genecounts"

IFS=$'_'
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
group="${PWD}/group"

# create two arrays for loop generation of grouping files
sampletype=("Clone1" "Clone2" "WT")
sampletime=("0" "24" "48")
# treatment=("Uninduced" "Induced")
for item1 in "${sampletype[@]}"
do
    for item2 in "${sampletime[@]}"
    do
	cat ${PWD}/fastq/Tco2.fqfiles |grep -w "$item1" | grep -w "$item2" | cut -f1 >> ${group}/"$item1"_"$item2".temp
	cat ${group}/"$item1"_"$item2".temp | awk '{gsub("Tco", "Tco-"); print}' > ${group}/"$item1"_"$item2".txt
	rm -f ${group}/*temp
	done
done

# delete empty files
find ${group}/ -type f -empty -delete

# 5.2 create gene descriptions file

# clear the previous results
rm -fr ${PWD}/average
# create a folder for grouping
mkdir ${PWD}/average
average="${PWD}/average"

cut -f4,5 ${PWD}/TriTrypDB-46_TcongolenseIL3000_2019.bed > ${average}/gene_descriptions.txt

# 5.3 calculate the mean expression levels for each group and add gene descriptions

# 1) extract the sixth column of each bedtools output file
# 2) paste them together according to grouping outcomes
# 3) calculate the mean expression levels
# 4) add gene descriptions

for item in $(ls ${group})
do
touch ${average}/${item}.outputfile
	lines=()
	while read line
	do
	cut -f6 ${genecounts}/${line}.gene.counts.txt > ${average}/cut.${line}
	lines+=("$line")
	done < ${group}/${item}
	output="${average}/${item}.outputfile"
	for ((i = 0; i < ${#lines[@]}; i++))
	do
		currentfile="${average}/cut.${lines[i]}"
		paste "$currentfile" "$output" > ${average}/temp.txt
		mv ${average}/temp.txt ${output}
	done
	cat ${output} | source ${PWD}/average.sh > ${average}/avgtemp.${item}
	# add header for the
	echo -e "name\tdescriptions\taverage" > ${average}/avg.${item}
	paste ${average}/gene_descriptions.txt ${average}/avgtemp.${item} >>  ${average}/avg.${item}
done

rm -f ${average}/*temp* ${average}/*cut* ${average}/*output*

# 6 calculate the fold change -- compared to the WT control group over time

rm -fr ${PWD}/foldchange
# create a folder for grouping
mkdir ${PWD}/foldchange
foldchange="${PWD}/foldchange"

# paste the column together for further calculate
for item1 in "Clone1" "Clone2"
do 
	for item2 in "0" "24" "48"
	do
	paste <(cut -f3 ${average}/avg."${item1}"_"${item2}".txt) <(cut -f3 ${average}/avg.WT_"${item2}".txt) > ${foldchange}/"${item1}"-WT."${item2}"
	done
done

# calculate the foldchange compared to the wild type group
# sort the log2foldchange decreasing order
for item in $(ls $foldchange | grep "-")
do 
	cat ${foldchange}/${item} | source foldchange.sh > ${foldchange}/${item}.fc.temp
	paste ${average}/gene_descriptions.txt ${foldchange}/${item}.fc.temp > ${foldchange}/${item}.fc.txt
	sort -k3,3 -nr -t$'\t' ${foldchange}/${item}.fc.txt > ${foldchange}/${item}.fc.sorted
done

rm -f ${foldchange}/*temp
