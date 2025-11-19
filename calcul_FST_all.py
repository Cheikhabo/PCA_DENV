import allel
import numpy as np
import pandas as pd
import os
from itertools import combinations

VCF_FILE = "/home/abo/stage/dengue_virus/Illumina0/DENV/ALLDENV/diversite/fst_par_snp/merged_all_fixed.vcf.gz"
POP_DIR = "/home/abo/stage/dengue_virus/Illumina0/DENV/ALLDENV/diversite/fst_par_snp/populations_lists"
OUT_DIR = "/home/abo/stage/dengue_virus/Illumina0/DENV/ALLDENV/diversite/fst_par_snp/resultats"
FST_THRESHOLD = 0.25

os.makedirs(OUT_DIR, exist_ok=True)

# Lecture du VCF
print("Lecture du VCF...")
callset = allel.read_vcf(VCF_FILE, fields=['samples','calldata/GT','variants/CHROM','variants/POS','variants/REF','variants/ALT'])
samples = callset['samples']
gt = allel.GenotypeArray(callset['calldata/GT'])
chrom = callset['variants/CHROM']
pos = callset['variants/POS']
ref = callset['variants/REF']
alt = callset['variants/ALT']

# Construire indices populations
pop_files = [f for f in os.listdir(POP_DIR) if f.endswith('.txt')]
pop_idx = {}
for pop_file in pop_files:
    pop_name = os.path.splitext(pop_file)[0]
    with open(os.path.join(POP_DIR, pop_file)) as f:
        ids = [l.strip() for l in f if l.strip()]
    idx = []
    for s in ids:
        loc = np.where(samples == s)[0]
        if len(loc) == 1:
            idx.append(int(loc[0].item()))
    pop_idx[pop_name] = np.array(idx, dtype=int)

# Fonction pour allele counts haploïdes avec traitement NaN
def haploid_ac(popname):
    idx = pop_idx[popname]
    sub = gt[:, idx]
    # Convertir en n_alt haploïde
    ac = sub.to_n_alt()  # shape (n_variants, n_samples)
    # Remplacer les valeurs manquantes (-1) par 0
    ac = np.where(ac < 0, 0, ac)
    return ac

# Paires de populations
pop_names = list(pop_idx.keys())
pairs = list(combinations(pop_names, 2))

for pop1, pop2 in pairs:
    print(f"Calcul FST : {pop1} vs {pop2}")
    ac1 = haploid_ac(pop1)
    ac2 = haploid_ac(pop2)

    # Initialiser tableau FST
    fst_array = np.full(ac1.shape[0], np.nan)  # taille n_variants

    # Calcul SNP par SNP
    for i in range(ac1.shape[0]):
        a1 = ac1[i, :].reshape(-1,1)
        a2 = ac2[i, :].reshape(-1,1)
        try:
            fst, _, _ = allel.hudson_fst(a1, a2)
            fst_array[i] = fst[0]  # hudson_fst retourne array
        except Exception:
            fst_array[i] = np.nan  # si impossible de calculer

    # DataFrame
    df = pd.DataFrame({
        'CHR': chrom,
        'POS': pos,
        'REF': ref,
        'ALT': [','.join(a) if isinstance(a, (list,tuple)) else str(a) for a in alt],
        'FST': fst_array
    })

    # Fichier complet
    out_file = os.path.join(OUT_DIR, f"{pop1}_vs_{pop2}_fst.csv")
    df.to_csv(out_file, index=False)

    # SNP high-FST
    df_high = df[df['FST'] >= FST_THRESHOLD]
    out_file_high = os.path.join(OUT_DIR, f"{pop1}_vs_{pop2}_highFST.csv")
    df_high.to_csv(out_file_high, index=False)

    print(f"Fichiers générés : {out_file} et {out_file_high}")

print("Tous les calculs FST haploïdes par SNP sont terminés !")
