set -e

bamfile=$1 # bamfile, without the final extension (.bam)
platform=$2 # what was the sequencing platform used? For RG addition
library=$3 # what was the sequencing library? For RG addition
sampleID=$4 # what is the sample ID? For RG addition
flowcell=$5 # what is the flowcell ID?
lane=$6 # what is the lane?

bn=$(basename "$bamfile" .bam)
dir=$(dirname "$bamfile")

picard AddOrReplaceReadGroups \
        -I $bamfile \
        -O ${dir}/${bn}.RG.bam \
        -RGLB ${library} \
        -RGPL ${platform} \
        -RGPU ${flowcell}.${lane} \
        -RGSM ${sampleID} \
        -VALIDATION_STRINGENCY LENIENT

samtools index ${dir}/${bn}.RG.bam
