#!/bin/bash


FASTQ_DIR="data"
OUTPUT_DIR="CHIKV_out"    
REFERENCE_GENOME="reference_Chikungunya.fasta"

mkdir -p "$OUTPUT_DIR/fastqc_results"
mkdir -p "$OUTPUT_DIR/trimmed"
mkdir -p "$OUTPUT_DIR/alignment"
mkdir -p "$OUTPUT_DIR/quast"
mkdir -p "$OUTPUT_DIR/variants"
mkdir -p "$OUTPUT_DIR/consensus"

#echo "Running FastQC on all FASTQ files..."
#fastqc $FASTQ_DIR/*.fastq.gz -o $OUTPUT_DIR/fastqc_results

echo "Trimming reads with Trimmomatic..."
for sample in data/*_L001_R1_001.fastq.gz; do
  sample=$(basename $sample)
  sample_name="${sample%_L*}"

  R1="$FASTQ_DIR/${sample_name}_L001_R1_001.fastq.gz"
  R2="$FASTQ_DIR/${sample_name}_L001_R2_001.fastq.gz"


  R1_trimmed="$OUTPUT_DIR/trimmed/${sample_name}_R1_trimmed.fastq.gz"
  R2_trimmed="$OUTPUT_DIR/trimmed/${sample_name}_R2_trimmed.fastq.gz"
 
  java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 "$R1" "$R2" "$R1_trimmed" "$OUTPUT_DIR/trimmed/${sample_name}_R1_unpaired.fastq.gz" \
                "$R2_trimmed" "$OUTPUT_DIR/trimmed/${sample_name}_R2_unpaired.fastq.gz" \
                 ILLUMINACLIP:/home/biobacteria/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36


  echo "Aligning reads with BWA..."
    
  SAM_FILE="$OUTPUT_DIR/alignment/${sample_name}.sam"

  bwa mem -t 12 "$REFERENCE_GENOME" "$R1_trimmed" "$R2_trimmed" > "$SAM_FILE"

  echo "Processing SAM files with SAMtools..."
  BAM_FILE="$OUTPUT_DIR/alignment/${sample_name}.bam"
  SORTED_BAM_FILE="$OUTPUT_DIR/alignment/${sample_name}_sorted.bam"
  
  samtools view -bS "$SAM_FILE" > "$BAM_FILE"
  
  samtools sort "$BAM_FILE" -o "$SORTED_BAM_FILE"
  
  samtools index "$SORTED_BAM_FILE"
  
#  echo "Calling variants with FreeBayes..."

  VCF_FILE="$OUTPUT_DIR/variants/${sample_name}_variants.vcf"
  VCF_FILE_GZ="$OUTPUT_DIR/variants/${sample_name}_variants.vcf.gz"
  
#  freebayes -f "$REFERENCE_GENOME" "$SORTED_BAM_FILE" > "$VCF_FILE"

  echo "Generating consensus sequences..."
  
  CONSENSUS_FILE="$OUTPUT_DIR/consensus/${sample_name}_consensus.fasta"
  
  bcftools mpileup -f "$REFERENCE_GENOME" "$SORTED_BAM_FILE" | bcftools call -mv -Oz -o "$VCF_FILE"
  bgzip -c "$VCF_FILE" > "$VCF_FILE_GZ"
  bcftools index "$VCF_FILE_GZ"
  bcftools consensus -f "$REFERENCE_GENOME" "$VCF_FILE_GZ" > "$CONSENSUS_FILE"
  
  echo "Running QUAST on mapped reads for $sample_name..."
  quast $CONSENSUS_FILE --bam "$SORTED_BAM_FILE" -r "$REFERENCE_GENOME" -o "$OUTPUT_DIR/quast/${sample_name}" 


done

echo "Pipeline completed successfully!"
