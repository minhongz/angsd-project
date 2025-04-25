#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=download
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=62G
#SBATCH --output="./logs/download.out"
#SBATCH --error="./logs/download.err"

sample_id=(SRR6244046 SRR6244047 SRR6244048 SRR6244049 SRR6244050 SRR6244051 \
           SRR6244052 SRR6244053 SRR6244054 SRR6244055 SRR6244056 SRR6244057 SRR6244058)

sample_name=(bm_cd69neg_ind1_1 bm_cd69neg_ind1_2 bm_cd69neg_ind2 bm_cd69neg_ind3 ln_cd69neg_ind1 \
             ln_cd69neg_ind3 ln_cd69neg_ind2 bm_cd69pos_ind1 bm_cd69pos_ind2 bm_cd69pos_ind3 \
             ln_cd69pos_ind1 ln_cd69pos_ind3 ln_cd69pos_ind2)

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
done

# === Merge bm_cd69neg_ind1 ===
sample="bm_cd69neg_ind1"
fq1="./raw_fastq/${sample}_1.fastq.gz"
fq2="./raw_fastq/${sample}_2.fastq.gz"

split_samples=("bm_cd69neg_ind1_1" "bm_cd69neg_ind1_2")
s1_fq1="./raw_fastq/${split_samples[0]}_1.fastq.gz"
s1_fq2="./raw_fastq/${split_samples[0]}_2.fastq.gz"
s2_fq1="./raw_fastq/${split_samples[1]}_1.fastq.gz"
s2_fq2="./raw_fastq/${split_samples[1]}_2.fastq.gz"

echo "Merging:"
echo "$s1_fq1 + $s2_fq1 -> $fq1"
echo "$s1_fq2 + $s2_fq2 -> $fq2"

zcat $s1_fq1 $s2_fq1 | gzip > $fq1
zcat $s1_fq2 $s2_fq2 | gzip > $fq2
fastqc "$fq1" "$fq2" -o ./raw_fastq/fastqc_raw_report