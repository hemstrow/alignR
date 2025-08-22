set -e

vcf=$1
ref=$2
QD="${3}"
FS="${4}"
SOR="${5}"
MQ="${6}"
MQRS="${7}"
RPRS="${8}"
GQ="${9}"
mem=${10}


bn=$(basename "$vcf" .gz)
bn=$(basename "$bn" .vcf)
dir=$(dirname "$vcf")

gatk --java-options "-Xmx${mem}g -Xms${mem}g" VariantFiltration \
       -R $ref \
       -V $vcf \
       --filter-name "QDf" \
       --filter-expression "QD < $QD" \
       --filter-name "FSf" \
       --filter-expression "FS > $FS" \
       --filter-name "SORf" \
       --filter-expression "SOR > $SOR" \
       --filter-name "MQf" \
       --filter-expression "MQ < $MQ" \
       --filter-name "MQRSf" \
       --filter-expression "MQRankSum < $MQRS" \
       --filter-name "RPRSf" \
       --filter-expression "ReadPosRankSum < $RPRS" \
       -O ${dir}/${bn}_hard_filt.vcf

bcftools view --types snps -m 2 -M 2 ${dir}/${bn}_hard_filt.vcf > ${dir}/${bn}_hard_filt_temp.vcf

vcftools --vcf ${dir}/${bn}_hard_filt_temp.vcf \
        --remove-filtered-all \
        --minGQ $GQ \
        --recode \
        --recode-INFO-all \
        --out ${dir}/${bn}_hard_filt_pass

rm ${dir}/${bn}_hard_filt_temp.vcf
