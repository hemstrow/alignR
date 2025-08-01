set -e

c1=$1 # A file containing the names of the files containing data for each individual. The .fastq at the end is cutt off.
ref=$2 # the path to the genome for the alignment

bwa mem $ref ${c1} | samtools view -Sb - | samtools sort -o ${c1}.sort.flt.bam -

samtools index ${c1}.sort.flt.bam
samtools flagstat ${c1}.sort.flt.bam.flagstat
samtools stats ${c1}.sort.flt.bam > ${c1}.sort.flt.bam.stats
