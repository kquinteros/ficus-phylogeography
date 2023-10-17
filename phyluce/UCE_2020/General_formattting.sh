#!/bin/bash

#Remove prefixes from Feb2020 fastq files in raw_fastq
for file in Lab-ID-*; do mv "$file" "${file#Lab-ID-}";done;

#Remove text from May2018 fastq files in raw_fastq
for file_name in FW*.gz
do
  new_file_name=$(sed 's/_[^_]*\_/_/g' <<< "$file_name");
  mv "$file_name" "$new_file_name";
done

#Split assembly_pet_Feb2020.conf into mulitple files
split -a 2 -l 6 assembly_pet_Feb2020.conf assembly_Feb2020_

for file in assembly_Feb2020_*
do
  mv "$file" "${file}.conf"
done

#create tmp files with header
vi tmp_header

for file in assembly_Feb2020_*
do
    head -n 2 tmp_header > tmp_file
    cat "$file" >> tmp_file
    mv -f tmp_file "$file"
done
rm tmp_header

#Create cmd files
ls assembly_conf > tmp.txt
#edit in text editor and replace assembly_Feb2020 with cmd_trinity
# and .conf with .#!/bin/sh
input="tmp.txt"
while IFS= read -r line
do
  touch "$line"
done < "$input"
rm tmp.txt

#create tmp file with sbatch commands
#copy slurm commands into each file
for i in cmd_trinity*.sh
do
   cat tmp > $i
done


for i in cmd_spades*.sh
do
   cat tmp > $i
done
