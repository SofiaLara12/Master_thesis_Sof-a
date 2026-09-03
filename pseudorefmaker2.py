import os
import sys
from pathlib import Path
from Bio import SeqIO

ARCHIVO_TSV = "captus-extract_stats.tsv"
ARCHIVO_DONANTES = "donantes_pseudoreferencia.txt"
OUTPUT_FASTA = "pseudoreferencia_taxonomica.fasta"
DIR_EXTRACTIONS = Path(".")

if not Path(ARCHIVO_DONANTES).exists():
    print(f"[!] Error: No existe '{ARCHIVO_DONANTES}'.")
    sys.exit(1)

muestra_principal = None
donantes_map = {}

with open(ARCHIVO_DONANTES, "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        if parts[0] == "MUESTRA_PRINCIPAL":
            muestra_principal = parts[1]
        elif parts[0] == "DONANTE":
            donante = parts[1]
            loci = parts[3].split(",") if len(parts) > 3 else []
            donantes_map[donante] = set(loci)

locus_a_muestra = {}
for donante, loci_set in donantes_map.items():
    for locus in loci_set:
        locus_a_muestra[locus] = donante

print(f"[+] Muestra principal: {muestra_principal}")
print(f"[+] Donantes asignados: {len(donantes_map)}")

mapa_contigs = {}

print(f"[+] Leyendo mapa de contigs desde {ARCHIVO_TSV}...")

with open(ARCHIVO_TSV, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if parts[0] == "sample_name":
            continue

        sample = parts[0]
        locus = parts[2]
        ctg_names = parts[22] if len(parts) > 22 else ""

        if ctg_names and ctg_names != "NA":
            primer_ctg = ctg_names.split(",")[0].strip()
            primer_ctg_id = primer_ctg.split()[0]
            key = (sample, locus)
            if key not in mapa_contigs:
                mapa_contigs[key] = primer_ctg_id

print(f"[+] Total de combinaciones Muestra-Locus mapeadas: {len(mapa_contigs)}")

target_contigs = {}


for (sample, locus), ctg_id in mapa_contigs.items():
    if sample == muestra_principal:
        target_contigs[ctg_id] = (locus, muestra_principal)

for donante, loci_set in donantes_map.items():
    for locus in loci_set:
        key = (donante, locus)
        if key in mapa_contigs:
            ctg_id = mapa_contigs[key]
            target_contigs[ctg_id] = (locus, donante)

print(f"[+] Total de secuencias objetivo a extraer: {len(target_contigs)}")

print("\n[+] Indexando y extrayendo de los archivos FASTA...")

fasta_files = list(DIR_EXTRACTIONS.rglob("*.fasta")) + \
              list(DIR_EXTRACTIONS.rglob("*.FNA")) + \
              list(DIR_EXTRACTIONS.rglob("*.fa"))

secuencias_extraidas = {}

for fasta_path in fasta_files:
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            record_id = record.id.split()[0]
            if record_id in target_contigs:
                locus, muestra = target_contigs[record_id]
                if locus not in secuencias_extraidas:
                    secuencias_extraidas[locus] = (muestra, str(record.seq))
    except Exception:
        continue

print(f"\n[+] Recolección finalizada. Secuencias recuperadas: {len(secuencias_extraidas)}")

with open(OUTPUT_FASTA, "w") as out:
    for locus_id, (sample, seq) in sorted(secuencias_extraidas.items()):
        out.write(f">{locus_id}|{sample}\n{seq}\n")

print(f"[✓] Archivo '{OUTPUT_FASTA}' generado con éxito.")
