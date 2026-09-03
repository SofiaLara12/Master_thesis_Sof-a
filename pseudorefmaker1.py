import sys
from collections import defaultdict

ARCHIVO = "captus-extract_stats.tsv"
MUESTRA_OBJETIVO = "Microcycas"

muestras_loci = defaultdict(set)
todos_los_loci = set()

with open(ARCHIVO, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        parts = line.strip().split("\t")

        if parts[0] == "sample_name":
            continue

        sample = parts[0]
        locus = parts[2]

        muestras_loci[sample].add(locus)
        todos_los_loci.add(locus)

target_sample = None
for s in muestras_loci.keys():
    if MUESTRA_OBJETIVO.lower() in s.lower():
        target_sample = s
        break

if not target_sample:
    print(f"[!] No se encontró ninguna muestra que contenga '{MUESTRA_OBJETIVO}'")
    sys.exit(1)

loci_microcycas = muestras_loci[target_sample]
loci_faltantes = todos_los_loci - loci_microcycas

print(f"=== RESULTADOS (IDÉNTICOS A R) ===")
print(f"Muestra principal: {target_sample}")
print(f"Total loci en el dataset: {len(todos_los_loci)}")
print(f"Loci recuperados por Microcycas: {len(loci_microcycas)}")
print(f"Loci faltantes en Microcycas:   {len(loci_faltantes)}\n")

loci_por_cubrir = set(loci_faltantes)
mapa_donantes = {
    s: loci for s, loci in muestras_loci.items() if s != target_sample
}

donantes_seleccionados = []

while loci_por_cubrir:
    mejor_donante = None
    max_aporte = set()

    for muestra, loci in mapa_donantes.items():
        aporte = loci.intersection(loci_por_cubrir)
        if len(aporte) > len(max_aporte):
            max_aporte = aporte
            mejor_donante = muestra

    if not mejor_donante or len(max_aporte) == 0:
        print(
            f"[!] Advertencia: {len(loci_por_cubrir)} loci no están presentes en NINGUNA otra muestra."
        )
        break

    donantes_seleccionados.append(
        (mejor_donante, len(max_aporte), sorted(list(max_aporte)))
    )
    loci_por_cubrir -= max_aporte
    del mapa_donantes[mejor_donante]

print(f"=== SELECCIÓN ÓPTIMA DE MUESTRAS DONANTES ===")
print(
    f"Se necesitan {len(donantes_seleccionados)} muestra(s) adicional(es):\n"
)

for i, (donante, num_aportados, lista_loci) in enumerate(
    donantes_seleccionados, 1
):
    print(f"  {i}. {donante}")
    print(f"     Aporta: {num_aportados} loci faltantes")
    print(
        f"     Ejemplos de loci: {lista_loci[:3]}{'...' if len(lista_loci) > 3 else ''}\n"
    )


with open("donantes_pseudoreferencia.txt", "w") as f:
    f.write(f"MUESTRA_PRINCIPAL\t{target_sample}\t{len(loci_microcycas)}\n")
    for donante, num, loci_list in donantes_seleccionados:
        f.write(f"DONANTE\t{donante}\t{num}\t{','.join(loci_list)}\n")

print(
    "[+] Archivo 'donantes_pseudoreferencia.txt' generado con éxito para el siguiente paso."
)
