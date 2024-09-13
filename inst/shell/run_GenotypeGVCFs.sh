set -e

bedfile=$1
ref=$2
mem=$3
tmp_dir=$4
javapath=$5
gatkpath=$6

bn=$(basename "$bedfile" .bed)
dir=$(dirname "$bedfile")

# run the genotyping
$javapath -jar -Xmx${mem}g -Xms${mem}g $gatkpath GenotypeGVCFs \
  -R $ref \
  -L $bedfile \
  -V gendb://${dir}/${bn}_db \
  -O ${dir}/${bn}_raw.vcf \
  --tmp-dir $tmp_dir \
  --new-qual \
  --max-alternate-alleles 2 \
  --disable-bam-index-caching
