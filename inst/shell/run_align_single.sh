set -e

c1=$1 # A file containing the names of the files containing data for each individual. The .fastq at the end is cutt off.
ref=$3 # the path to the genome for the alignment

bwa mem $ref ${c1} | samtools view -Sb - | samtools sort -o ${c2}.sort.bam -
samtools index ${c2}.sort.flt.bam
