set -e

conda=$1
name=$2

# not running things sequentially can cause big errors on SOME machines, like my WSL install.
$conda init bash
$conda create -n $name -y
$conda activate $name
$conda install -c bioconda -c conda-forge -y r-base
$conda install -c bioconda -c conda-forge -y samtools
$conda install -c bioconda -c conda-forge -y bamtools
$conda install -c bioconda -c conda-forge -y vcftools
$conda install -c bioconda -c conda-forge -y bcftools
$conda install -c bioconda -c conda-forge -y bwa
$conda install -c bioconda -c conda-forge -y perl
$conda install -c bioconda -c conda-forge -y stacks
$conda install -c bioconda -c conda-forge -y gatk4
$conda install -c bioconda -c conda-forge -y picard

R -e "install.packages('remotes')"
R -e "remotes::install_github('hemstrow/alignR')"

$conda deactivate
