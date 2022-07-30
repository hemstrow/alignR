list=${1}
minInd=$2
genotyper=$3
SNP_pval=$4
doGeno=$5
postCutoff=$6
minQ=$7
minMapQ=$8
outfile=$9
doVcf=${10}
par=${11}

angsd -bam ${list} \
  -GL $genotyper \
  -out $outfile \
  -doMaf 2 \
  -doMajorMinor 1 \
  -SNP_pval $SNP_pval \
  -doGeno $doGeno \
  -doPost 2 \
  -postCutoff $postCutoff \
  -minQ $minQ \
  -minMapQ $minMapQ \
  -minInd $minInd \
  -doBcf $doVcf \
  -nThreads $par

