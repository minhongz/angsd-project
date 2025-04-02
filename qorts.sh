#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=qorts_all
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=62G
#SBATCH --output="qorts_all.out"
#SBATCH --error="qorts_all.err"

conda activate qorts
qorts_jar="/athena/angsd/scratch/diz4009/qorts/hartleys-QoRTs-099881f/QoRTs.jar"

sample_name=(bm_cd69neg_ind1 bm_cd69neg_ind2 bm_cd69neg_ind3 ln_cd69neg_ind1 \
             ln_cd69neg_ind3 ln_cd69neg_ind2 bm_cd69pos_ind1 bm_cd69pos_ind2 \
             bm_cd69pos_ind3 ln_cd69pos_ind1 ln_cd69pos_ind3 ln_cd69pos_ind2)

for i in "${!sample_name[@]}"; do
    sample="${sample_name[$i]}"

    java -Xmx60G -jar "$qorts_jar" QC \
        --generatePdfReport \
        --stranded \
        --numThreads 24 \
        ./alignment/${sample}/*.bam \
        ./refGenomes/gencode.v47.basic.annotation.gtf \
        ./qorts/${sample}/
done

# sample=bm_cd69pos_ind3
# mkdir -p ./qorts
# java -Xmx60G -jar "$qorts_jar" QC \
#     --generatePdfReport \
#     --stranded \
#     --numThreads 24 \
#     ./alignment/${sample}/*.bam \
#     ./refGenomes/gencode.v47.basic.annotation.gtf \
#     ./qorts/${sample}/