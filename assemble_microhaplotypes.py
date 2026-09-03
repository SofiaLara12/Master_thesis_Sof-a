#!/usr/bin/env python3

import csv
import os
import argparse
from collections import defaultdict

def combine_microhaplotypes_from_csv(input_csv, output_file):
    """Combines phase-contiguous microhaplotypes from a CSV file and outputs the result."""

    if not os.path.exists(input_csv):
        print(f"ERROR: Input CSV file not found: {input_csv}")
        return

    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    haplotypes = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {
        "hap1": "",
        "hap2": "",
        "variant_positions": []
    })))

    with open(input_csv, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene_name = row["Gene"]
            sample_id = row["Sample"]
            block_id = row["Haplotype_Block"]


            if block_id == "Unknown":
                continue
            if row.get("Reason_Unknown", "").strip():
                continue

            genotype = row["Genotype"].split("|")
            if len(genotype) != 2:
                continue
            allele1 = row["Allele1"]
            allele2 = row["Allele2"]

            try:
                variant_pos = int(row["Variant_Pos"])
            except ValueError:
                continue

           haplotype_data = haplotypes[gene_name][sample_id][block_id]

            haplotype_data["hap1"] += allele1
            haplotype_data["hap2"] += allele2
            haplotype_data["variant_positions"].append(variant_pos)

    with open(output_file, "w", newline="") as csvfile:
        fieldnames = [
            "Gene", "Sample", "Haplotype_Block", 
            "First_SNP_Position", "Last_SNP_Position", "Haplotype_Block_Length",
            "Microhaplotype_1", "Microhaplotype_2", "Microhaplotype_Length"
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for gene in haplotypes:
            for sample in haplotypes[gene]:
                for block in haplotypes[gene][sample]:
                    hap1 = haplotypes[gene][sample][block]["hap1"]
                    hap2 = haplotypes[gene][sample][block]["hap2"]
                    variant_positions = haplotypes[gene][sample][block]["variant_positions"]

                    microhap_length = (len(hap1) + len(hap2)) / 2 if hap1 and hap2 else 0

                    if variant_positions:
                        first_snp_pos = min(variant_positions)
                        last_snp_pos = max(variant_positions)
                        haplotype_block_length = last_snp_pos - first_snp_pos if last_snp_pos != first_snp_pos else 1
                    else:
                        first_snp_pos = last_snp_pos = haplotype_block_length = None

                    if round(microhap_length, 2) == 1:
                        continue

                    writer.writerow({
                        "Gene": gene,
                        "Sample": sample,
                        "Haplotype_Block": block,
                        "First_SNP_Position": first_snp_pos,
                        "Last_SNP_Position": last_snp_pos,
                        "Haplotype_Block_Length": round(haplotype_block_length, 2),
                        "Microhaplotype_1": hap1,
                        "Microhaplotype_2": hap2,
                        "Microhaplotype_Length": round(microhap_length, 2)
                    })

    print(f"Microhaplotype combination complete. Output saved to: {output_file}")
    print("icrohaplotypes assembled successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine phase-contiguous microhaplotypes from CSV.")
    parser.add_argument("-i", "--input_csv", required=True, help="Path to the input CSV file with phased genotypes.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output CSV file for combined microhaplotypes.")

    args = parser.parse_args()
    combine_microhaplotypes_from_csv(args.input_csv, args.output_file)
