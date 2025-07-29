set -e

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
doRf=${11}
rf=${12}
par=${13}
minMaf=${14}

if [[ $doRf -eq 1 ]]
then
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
  -minMaf $minMaf \
  -doBcf $doVcf \
  -rf $rf \
  -nThreads $par
else
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
  -minMaf $minMaf \
  -doBcf $doVcf \
  -nThreads $par
fi



