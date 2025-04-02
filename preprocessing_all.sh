#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=preprocessing_all
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=62G
#SBATCH --output="preprocessing_all.out"
#SBATCH --error="preprocessing_all.err"

sample_id=(SRR6244046 SRR6244047 SRR6244048 SRR6244049 SRR6244050 SRR6244051 \
           SRR6244052 SRR6244053 SRR6244054 SRR6244055 SRR6244056 SRR6244057 SRR6244058)

sample_name=(bm_cd69neg_ind1_1 bm_cd69neg_ind1_2 bm_cd69neg_ind2 bm_cd69neg_ind3 ln_cd69neg_ind1 \
             ln_cd69neg_ind3 ln_cd69neg_ind2 bm_cd69pos_ind1 bm_cd69pos_ind2 bm_cd69pos_ind3 \
             ln_cd69pos_ind1 ln_cd69pos_ind3 ln_cd69pos_ind2)

skip_samples=("bm_cd69neg_ind1_1" "bm_cd69neg_ind1_2")

for i in "${!sample_id[@]}"; do
    srr="${sample_id[$i]}"
    sample="${sample_name[$i]}"
    fq1="./raw_fastq/${srr}_1.fastq.gz"
    fq2="./raw_fastq/${srr}_2.fastq.gz"
    renamed_fq1="./raw_fastq/${sample}_1.fastq.gz"
    renamed_fq2="./raw_fastq/${sample}_2.fastq.gz"

    echo "Processing $srr -> $sample"

    # === Check if FASTQ files already downloaded ===
    conda activate angsd
    if [[ -f "$fq1" && -f "$fq2" ]]; then
        echo "FASTQ files for $srr already downloaded."

        # Rename if not already renamed
        if [[ ! -f "$renamed_fq1" ]]; then
            mv "$fq1" "$renamed_fq1"
            echo "Renamed $fq1 to $renamed_fq1"
        fi
        if [[ ! -f "$renamed_fq2" ]]; then
            mv "$fq2" "$renamed_fq2"
            echo "Renamed $fq2 to $renamed_fq2"
        fi
    elif [[ -f "$renamed_fq1" && -f "$renamed_fq2" ]]; then
        echo "Renamed FASTQ files for $sample already exist. Skipping download."
    else
        echo "Downloading FASTQ files for $srr"
        conda activate angsd
        # fastq-dump --split-files --gzip -O ./raw_fastq "$srr"
        fasterq-dump $srr --split-files -O ./raw_fastq --threads 16

        # Compress
        pigz -p 16 ./raw_fastq/${srr}_1.fastq
        pigz -p 16 ./raw_fastq/${srr}_2.fastq

        # Rename after download
        mv "$fq1" "$renamed_fq1"
        mv "$fq2" "$renamed_fq2"

        # === FastQC ===
        fastqc "$renamed_fq1" "$renamed_fq2" -o ./raw_fastq/fastqc_raw_report

    fi

    # === Skip samples from bw_cd69neg_ind1 ===
    if [[ " ${skip_samples[@]} " =~ " ${sample} " ]]; then
        echo "Skipping alignment for $sample (marked for later processing)"
        continue
    fi

    # === STAR Alignment ===
    mkdir -p "./alignment/${sample}"
    STAR --runMode alignReads \
        --runThreadN 16 \
        --genomeDir /athena/angsd/scratch/miz4010/project/refGenomes/hg38_STARindex \
        --readFilesIn "$renamed_fq1" "$renamed_fq2" \
        --readFilesCommand zcat \
        --alignIntronMax 100000 \
        --alignSJDBoverhangMin 3 \
        --outFileNamePrefix "./alignment/${sample}/${sample}.STAR." \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes All

    samtools index -@ 16 "./alignment/${sample}/${sample}.STAR.Aligned.sortedByCoord.out.bam"

    # === Qualimap ===
    conda activate qualimap
    qualimap bamqc -bam "./alignment/${sample}/${sample}.STAR.Aligned.sortedByCoord.out.bam" \
        -outformat PDF -outfile "${sample}.STAR.report.pdf"

    echo "Finished processing $sample"
done