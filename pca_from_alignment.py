from Bio import AlignIO
import numpy as np
from sklearn.manifold import MDS
import pandas as pd
import sys

# Nom fichier aligné (multi-FASTA)
input_fasta = "aligned_consensus.fasta"

print(f"Chargement de l'alignement depuis {input_fasta}...")
aln = AlignIO.read(input_fasta, "fasta")

n = len(aln)
names = [rec.id for rec in aln]

print(f"{n} séquences chargées, longueur alignement = {aln.get_alignment_length()}")

# Calcul matrice de distances p-distance
print("Calcul de la matrice de distances p-distance...")
D = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        diffs = 0
        valid_sites = 0
        for a, b in zip(aln[i].seq, aln[j].seq):
            if a == "-" or b == "-":
                continue
            if a.upper() == "N" or b.upper() == "N":
                continue
            valid_sites += 1
            if a != b:
                diffs += 1
        dist = diffs / valid_sites if valid_sites > 0 else 0
        D[i, j] = D[j, i] = dist

print("Matrice de distances calculée.")

# PCA via Classical MDS avec n_init fixé à 4 pour éviter le warning
print("Calcul du PCA via MDS...")
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42, n_init=4)
coords = mds.fit_transform(D)

# Sauvegarde résultats
df = pd.DataFrame(coords, index=names, columns=["PC1", "PC2"])
out_csv = "PCA_coords_from_distance.csv"
df.to_csv(out_csv)
print(f"Coordonnées PCA sauvegardées dans {out_csv}")
