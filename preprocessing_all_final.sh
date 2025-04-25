#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=preprocessing_all_final
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=62G
#SBATCH --output="./logs/preprocessing_all_final.out"
#SBATCH --error="./logs/preprocessing_all_final.err"

sample_name=(bm_cd69neg_ind2 bm_cd69neg_ind3 ln_cd69neg_ind1 \
             ln_cd69neg_ind2 ln_cd69neg_ind3 bm_cd69pos_ind1 bm_cd69pos_ind2 \
             bm_cd69pos_ind3 ln_cd69pos_ind2 ln_cd69pos_ind3)

for sample in "${sample_name[@]}"; do
    
    raw_fq_1="./raw_fastq/${sample}_1.fastq.gz"
    raw_fq_2="./raw_fastq/${sample}_2.fastq.gz"
    trimmed_fq_1="./trimmed_fastq/${sample}_1_trimmed.fastq.gz"
    trimmed_fq_2="./trimmed_fastq/${sample}_2_trimmed.fastq.gz"

    # === Remove rRNA ===
    conda activate angsd

    mkdir -p "./alignment_rRNA/${sample}"
    STAR --runMode alignReads \
        --genomeDir ./refGenomes/rRNA_STARindex \
        --readFilesIn $raw_fq_1 $raw_fq_2 \
        --readFilesCommand zcat \
        --alignIntronMax 1 \
        --runThreadN 24 \
        --outFileNamePrefix "./alignment_rRNA/${sample}/${sample}.rRNA." \
        --outSAMtype None \
        --outReadsUnmapped Fastx

    mkdir -p "./trimmed_fastq"
    mv ./alignment_rRNA/${sample}/${sample}.rRNA.Unmapped.out.mate1 ./trimmed_fastq/${sample}_1.trimmed.fastq
    mv ./alignment_rRNA/${sample}/${sample}.rRNA.Unmapped.out.mate2 ./trimmed_fastq/${sample}_2.trimmed.fastq
    pigz -p 24 ./trimmed_fastq/${sample}_1_trimmed.fastq
    pigz -p 24 ./trimmed_fastq/${sample}_2_trimmed.fastq
    
    mkdir -p ./trimmed_fastq/fastqc_trimmed_report
    fastqc "$trimmed_fq_1" "$trimmed_fq_2" -o ./trimmed_fastq/fastqc_trimmed_report

    # === STAR Alignment ===
    mkdir -p "./alignment_trimmed/${sample}"
    STAR --runMode alignReads \
        --runThreadN 24 \
        --genomeDir /athena/angsd/scratch/miz4010/project/refGenomes/hg38_STARindex \
        --readFilesIn "$trimmed_fq_1" "$trimmed_fq_2" \
        --readFilesCommand zcat \
        --alignIntronMax 100000 \
        --alignSJDBoverhangMin 3 \
        --outFileNamePrefix "./alignment_trimmed/${sample}/${sample}.STAR." \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes All

    samtools index -@ 24 "./alignment_trimmed/${sample}/${sample}.STAR.Aligned.sortedByCoord.out.bam"

    # === Qualimap ===
    conda activate qualimap
    qualimap bamqc -bam "./alignment_trimmed/${sample}/${sample}.STAR.Aligned.sortedByCoord.out.bam" \
        -outformat PDF -outfile "${sample}.STAR.report.pdf"

    # === Qorts ===
    conda activate qorts
    qorts_jar="/athena/angsd/scratch/diz4009/qorts/hartleys-QoRTs-099881f/QoRTs.jar"

    mkdir -p ./qorts_trimmed/${sample}/
    java -Xmx60G -jar "$qorts_jar" QC \
        --generatePdfReport \
        --stranded \
        --numThreads 24 \
        ./alignment_trimmed/${sample}/*.bam \
        ./refGenomes/gencode.v47.basic.annotation.gtf \
        ./qorts_trimmed/${sample}/

    echo "Finished processing $sample"
done

# === MultiQC ===
conda activate multiqc
multiqc ./alignment ./qorts_raw ./raw_fastq -o ./multiqc_reports/all_raw
multiqc ./alignment_trimmed ./qorts_trimmed ./trimmed_fastq -o ./multiqc_reports/all_trimmed

# === featureCounts ===
conda activate angsd
mkdir -p ./featureCounts
annot_file="./refGenomes/gencode.v47.basic.annotation.gtf"

featureCounts -T 24 \
    -a ${annot_file} \
    -o ./featureCounts/fc_output_trimmed_exon.txt \
    -s 2 \
    -t exon \
    -g gene_name \
    -p --countReadPairs \
    ./alignment_trimmed/*_cd69*/*.bam
    
