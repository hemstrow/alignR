set -e

bamlist=$1
ref=$2
genotyper=$3
SNP_pval=$4
minMapQ=$5
minQ=$6
par=$7

### Calculate paralog probabilities and get a list of paralogous loci  #
# locate snps to check
echo "Beginning SNP location for ${bamlist}.\n"
angsd -bam $bamlist -out results_snps_${bamlist} -ref $ref -GL $genotyper -doMajorMinor 1 -doMaf 2 -SNP_pval $SNP_pval -minMapQ $minMapQ -minQ $minQ -nThreads $par
gunzip results_snps_${bamlist}*.gz
cut -d$'\t' -f1-2  results_snps_${bamlist}.mafs | sed 1d > results_snps_${bamlist}.pos

# run mpileup and ngsParalog
echo "Beginning samtools mpileup for ${bamlist}.\n"
samtools mpileup -b $bamlist -l results_snps_${bamlist}.pos -f ${ref} > results_depth_${bamlist}
echo "Beginning ngsParalog calcLR for ${bamlist}.\n"
ngsParalog calcLR -infile results_depth_${bamlist} > results_paralogs_${bamlist}
