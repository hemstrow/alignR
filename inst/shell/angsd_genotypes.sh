list=${1}
minInd=$2
genotyper=$3
SNP_pval=$4
doGeno=$5
postCutoff=$6
minQ=$7
minMapQ=$8
outfile=$9
par=${10}

echo $outfile
echo $par

angsd -bam ${list} -GL $genotyper -out $outfile -doMaf 2 -doMajorMinor 1 -SNP_pval $SNP_pval -doGeno $doGeno -doPost 2 -postCutoff $postCutoff -minQ $minQ -minMapQ $minMapQ -minInd $minInd -nThreads $par
