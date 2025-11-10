import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Chemins vers les fichiers
pca_path = '../consensus_results/PCA_coords_from_distance.csv'
meta_path = '../metadone_finale_corrected.tsv'

# Chargement des données
pca = pd.read_csv(pca_path)
meta = pd.read_csv(meta_path, sep='\t')

# Affichage des colonnes (pour vérification)
print("Colonnes PCA:", pca.columns)
print("Colonnes Meta:", meta.columns)

# Fusion sur la colonne 'Echantillon'
df = pd.merge(pca, meta, on='Echantillon', how='left')

# Tracer le PCA en 2D coloré par location
plt.figure(figsize=(10, 8))
sns.scatterplot(data=df, x='PC1', y='PC2', hue='location', palette='tab10', s=100, edgecolor='black')

plt.title('')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend(title='Localite', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
