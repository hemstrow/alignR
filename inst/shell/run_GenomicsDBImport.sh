set -e

samplemap=$1
bedfile=$2
mem=$3
tmp_dir=$4
batchsize=$5

bn=$(basename "$bedfile" .bed)
dir=$(dirname "$bedfile")

gatk --java-options "-Xmx${mem}g -Xms${mem}g" GenomicsDBImport \
       --genomicsdb-workspace-path ${dir}/${bn}_db \
       --batch-size $batchsize \
       -L $bedfile \
       --sample-name-map $samplemap \
       --tmp-dir $tmp_dir

