#!/bin/bash
#we run this from the root folder where all the bams are located in subfolders
#it creates script
for fullname in `find ../no.filter/cancer_bed_files -name '*.bed'`
do
	name=$(basename "$fullname")
	echo $name
	./filter.bed.pl < $fullname > ../cancer_bed_files/$name
done
