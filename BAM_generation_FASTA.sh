#!/bin/bash
set -eo pipefail

#### Usage: bash BAM_generation_FASTA.sh pseudoreferencia.fasta SampleID

reference=$1
prefix=$2

fasta_input="${prefix}__captus-ext/01_coding_NUC/NUC_coding_NT.fna"

echo "=== Mapeando $fasta_input contra $reference ==="

bwa mem -t 4 "$reference" "$fasta_input" | samtools view -bS -q 20 -F 2048 -F 256 - | samtools sort - -o "${prefix}.sorted.bam"

gatk AddOrReplaceReadGroups \
  -I "${prefix}.sorted.bam" \
  -O "${prefix}.marked.noDups.bam" \
  -RGID "$prefix" \
  -RGLB lib1 \
  -RGPL illumina \
  -RGPU unit1 \
  -RGSM "$prefix"

samtools index "${prefix}.marked.noDups.bam"

rm "${prefix}.sorted.bam"


