#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
from collections import defaultdict
import pysam


def parse_args():
  parser = argparse.ArgumentParser(
      description=(
          "Extract microhaplotype SNPs from phased VCF for assembled"
          " transcriptomes."
      )
  )
  parser.add_argument(
      "-v",
      "--vcf",
      required=True,
      help="Input phased VCF file (uncompressed or bgzipped)",
  )
  parser.add_argument(
      "-g",
      "--gtf",
      required=True,
      help="GTF file containing block_id annotations",
  )
  parser.add_argument(
      "-i",
      "--imputation_log",
      required=True,
      help="CSV log of imputed genotypes",
  )
  parser.add_argument(
      "-o", "--output", required=True, help="Output CSV file"
  )
  return parser.parse_args()


def bgzip_and_index(vcf_path):
  """Compress and index a VCF if not already compressed."""
  if not vcf_path.endswith(".gz"):
    print(f"Compressing and indexing {vcf_path} ...")
    compressed_vcf = vcf_path + ".gz"
    with open(compressed_vcf, "wb") as f_out:
      subprocess.run(["bgzip", "-c", vcf_path], stdout=f_out, check=True)
    subprocess.run(["tabix", "-p", "vcf", compressed_vcf], check=True)
    return compressed_vcf
  else:
    index_path = vcf_path + ".tbi"
    if not os.path.exists(index_path):
      print(f"Index not found. Indexing {vcf_path} ...")
      subprocess.run(["tabix", "-p", "vcf", vcf_path], check=True)
    return vcf_path


def load_blocks_from_gtf(gtf_file):
  """Map positions from GTF blocks to block_ids."""
  blocks = defaultdict(list)
  with open(gtf_file) as f:
    for line in f:
      if line.startswith("#"):
        continue
      fields = line.strip().split("\t")
      chrom, _, feature, start, end, _, _, _, attributes = fields
      if feature != "block":
        continue
      start, end = int(start), int(end)
      block_id = attributes.split('block_id "')[1].split('"')[0]
      for pos in range(start, end + 1):
        blocks[chrom].append((pos, block_id))

  pos_to_block = defaultdict(str)
  for chrom in blocks:
    for pos, block_id in blocks[chrom]:
      pos_to_block[(chrom, pos)] = block_id
  return pos_to_block


def load_imputed_positions(log_file):
  """Load imputation log to determine which SNPs were imputed."""
  imputed = set()
  with open(log_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
      sample = row["Sample"]
      chrom_pos = row["SNP"]
      if ":" in chrom_pos:
        chrom, pos = chrom_pos.rsplit(":", 1)
        try:
          pos = int(pos)
        except ValueError:
          continue
        imputed.add((sample, chrom, pos))
  return imputed


def fetch_alleles(record, sample):
  """Extract alleles per individual."""
  gt = record.samples[sample].get("GT")
  if gt is None or "." in str(gt):
    return None, None, "./."
  alleles = [record.alleles[i] if i < len(record.alleles) else "N" for i in gt]
  genotype_str = "|".join(map(str, gt))
  return alleles[0], alleles[1], genotype_str


def main():
  args = parse_args()
  vcf_path = bgzip_and_index(args.vcf)
  vcf = pysam.VariantFile(vcf_path)

  pos_to_block = load_blocks_from_gtf(args.gtf)
  imputed = load_imputed_positions(args.imputation_log)

  print("VCF samples:", list(vcf.header.samples))
  print(f"Total imputed positions loaded: {len(imputed)}")

  with open(args.output, "w", newline="") as out_f:
    writer = csv.writer(out_f)
    header = [
        "Gene",
        "Variant_Pos",
        "Sample",
        "Genotype",
        "Allele1",
        "Allele2",
        "Haplotype_Block",
        "Imputed",
    ]
    writer.writerow(header)

    for record in vcf.fetch():
      chrom, pos = record.chrom, record.pos
      block_id = pos_to_block.get((chrom, pos), "NA")

      for sample in record.samples:
        allele1, allele2, genotype_str = fetch_alleles(record, sample)
        if allele1 is None:
          continue

        imputed_flag = "Yes" if (sample, chrom, pos) in imputed else "No"

        writer.writerow([
            chrom,
            pos,
            sample,
            genotype_str,
            allele1,
            allele2,
            block_id,
            imputed_flag,
        ])

  print(f"✅ Microhaplotype table exported to {args.output}")


if __name__ == "__main__":
  main()
