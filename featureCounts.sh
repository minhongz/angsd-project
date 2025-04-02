#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=featureCounts
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --output="featureCounts.out"
#SBATCH --error="featureCounts.err"

conda activate angsd

mkdir -p ./featureCounts
annot_file="./refGenomes/gencode.v47.basic.annotation.gtf"

featureCounts -T 24 \
    -a ${annot_file} \
    -o ./featureCounts/fc_output.txt \
    -s 2 \
    -t gene \
    -g gene_id \
    -M -O -p --countReadPairs \
    ./alignment/*_cd69*/*.bam

    