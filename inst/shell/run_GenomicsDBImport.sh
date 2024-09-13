set -e

samplemap=$1
bedfile=$2
mem=$3
tmp_dir=$4
javapath=$5
gatkpath=$6
batchsize=$7

bn=$(basename "$bedfile" .bed)
dir=$(dirname "$bedfile")

$javapath -jar -Xmx${mem}g -Xms${mem}g $gatkpath GenomicsDBImport \
       --genomicsdb-workspace-path ${dir}/${bn}_db \
       --batch-size $batchsize \
       -L $bedfile \
       --sample-name-map $samplemap \
       --tmp-dir $tmp_dir

