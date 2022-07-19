c1=$1 # A file containing the names of the files containing data for each individual (RA). The .fastq at the end is cutt off.
c2=$2 # A file containing the names of the files containing data for each individual (RA). The .fastq at the end is cutt off.
ref=$3 # the path to the genome for the alignment

bwa mem $ref ${c1}.fastq ${c2}.fastq | samtools view -Sb - | samtools sort - -n -o ${c1}.sort.bam # align and sort by name
samtools fixmate -r -m ${c1}.sort.bam ${c1}.fixmate.bam # fixmate
samtools sort -o ${c1}.psort.bam ${c1}.fixmate.bam # sort by position
samtools markdup -r ${c1}.psort.bam ${c1}.markdup.bam # remove dups
samtools view -q 5 -b ${c1}.markdup.bam > ${c1}.q1.bam # remove poorly mapped
samtools sort -n -o ${c1}.namesort.bam ${c1}.q1.bam # sort by name again
samtools fixmate -m ${c1}.namesort.bam ${c1}.fixmate.bam # filter bad mates again
samtools view -f 0x2 -b ${c1}.fixmate.bam > ${c1}.flt.bam # remove improper pairs
samtools sort -o ${c1}.sort.flt.bam ${c1}.flt.bam # sort
samtools index ${c1}.sort.flt.bam # index

# clean
rm ${c1}.markdup.bam
rm ${c1}.fixmate.bam
rm ${c1}.q1.bam
rm ${c1}.psort.bam
rm ${c1}.namesort.bam
