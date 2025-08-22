set -e

bedfile=$1
ref=$2
mem=$3
new_qual=$4
tmp_dir=$5

bn=$(basename "$bedfile" .bed)
dir=$(dirname "$bedfile")

# run the genotyping
if [[ $newqual -eq 1 ]]; then
  gatk --java-options "-Xmx${mem}g -Xms${mem}g" GenotypeGVCFs \
    -R $ref \
    -L $bedfile \
    -V gendb://${dir}/${bn}_db \
    -O ${dir}/${bn}_raw.vcf \
    --tmp-dir $tmp_dir \
    --new-qual \
    --max-alternate-alleles 2 \
    --disable-bam-index-caching
else
  gatk --java-options "-Xmx${mem}g -Xms${mem}g" GenotypeGVCFs \
    -R $ref \
    -L $bedfile \
    -V gendb://${dir}/${bn}_db \
    -O ${dir}/${bn}_raw.vcf \
    --tmp-dir $tmp_dir \
    --max-alternate-alleles 2 \
    --disable-bam-index-caching
fi


