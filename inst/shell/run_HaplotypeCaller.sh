set -e

bamfile=$1
reference=$2
tmp_dir=$3
mem=$4
min_base_quality_score=$5
minimum_mapping_quality=$6
region=$7


if [[ "$region" != "all" && "$region" != "" ]]; then
  if [[ -f "$region" ]]; then
    rprint=$(basename "$region" .list)
    rprint=$(basename "$rprint" .bed)
  else
    rprint="$region"
  fi

  gatk --java-options "-Xmx${mem}g -Xms${mem}g -Djava.io.tmpdir=${tmp_dir}" HaplotypeCaller \
        -ERC GVCF \
        -R $reference \
        -I $bamfile \
        --min-base-quality-score $min_base_quality_score \
        --minimum-mapping-quality $minimum_mapping_quality \
        -O ${bamfile}-${rprint}.hapcalls.gvcf.gz \
        -L $region

else
  gatk --java-options "-Xmx${mem}g -Xms${mem}g -Djava.io.tmpdir=${tmp_dir}" HaplotypeCaller \
          -ERC GVCF \
          -R $reference \
          -I $bamfile \
          --min-base-quality-score $min_base_quality_score \
          --minimum-mapping-quality $minimum_mapping_quality \
          -O ${bamfile}-all.hapcalls.gvcf.gz
fi


