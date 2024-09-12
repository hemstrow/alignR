ref=$1
bedfile=$2
tmp_dir=$3
mem=$4

basebed=$(basename ${bedfile})

echo "build database start at: `date`"
gatk --java-options "-Xmx${mem}g -Xms${mem}g" \
        GenotypeGVCFs \
  -R $ref \
  -L $bedfile \
  -V gendb://${basebed}_db \
  -O raw_${basebed}.vcf \
  --tmp-dir $tmp_dir


rm -r tmp_dir

echo "Genotypes done at: `date`"
echo "Filtering:"

### filter to ID high-quaity SNPs
gatk --java-options "-Xmx2g -Xms2g" \
       VariantFiltration \
       -R $ref \
       -V raw_${bedfile}.vcf \
       --filter-name "QDf" \
       --filter-expression "QD < 2.0" \
       --filter-name "FSf" \
       --filter-expression "FS > 60.0" \
       --filter-name "SORf" \
       --filter-expression "SOR > 3.0" \
       --filter-name "MQf" \
       --filter-expression "MQ < 40.0" \
       --filter-name "MQRSf" \
       --filter-expression "MQRankSum < -12.5" \
       --filter-name "RPRSf" \
       --filter-expression "ReadPosRankSum < -8.0" \
       -O hard_filt_${bedfile}.vcf

bcftools view --types snps -m 2 -M 2 hard_filt_${bedfile}.vcf > hard_filt_temp_${bedfile}.vcf

vcftools --vcf hard_filt_temp_${bedfile}.vcf \
        --remove-filtered-all \
        --remove-indels \
        --max-missing 0.2 \
        --recode \
        --recode-INFO-all \
        --out hard_filt_pass_${bedfile}

rm hard_filt_temp_${bedfile}.vcf
