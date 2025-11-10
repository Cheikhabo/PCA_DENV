#!/bin/bash

# Dossier contenant les VCF compressés
VCF_DIR="/home/abo/stage/dengue_virus/Illumina0/DENV/ALLDENV/vcfs"

# Dossier contenant les fichiers de référence FASTA
REF_DIR="/home/abo/stage/dengue_virus/Illumina0/DENV/reference"

# Dossier de sortie pour les consensus
OUT_DIR="/home/abo/stage/dengue_virus/Illumina0/DENV/ALLDENV/consensus_results"

# Créer le dossier de sortie s'il n'existe pas
mkdir -p "$OUT_DIR"

# Boucle sur tous les fichiers .vcf.gz (et pas .csi)
for vcf in "$VCF_DIR"/*.vcf.gz
do
  # Ignorer les fichiers d'index .csi
  [[ "$vcf" == *.csi ]] && continue

  base=$(basename "$vcf" .vcf.gz)

  # Détection du sérotype
  if [[ "$base" == *"DENV1"* ]]; then
    ref="$REF_DIR/DENV1.fasta"
  elif [[ "$base" == *"DENV2"* ]]; then
    ref="$REF_DIR/DENV2.fasta"
  elif [[ "$base" == *"DENV3"* ]]; then
    ref="$REF_DIR/DENV3.fasta"
  elif [[ "$base" == *"DENV4"* ]]; then
    ref="$REF_DIR/DENV4.fasta"
  else
    echo "⚠️  Sérotype inconnu pour $base — fichier ignoré."
    continue
  fi

  echo "➡️  Génération du consensus pour $base avec référence $(basename "$ref")"
  
  # Génération du consensus
  bcftools consensus -f "$ref" "$vcf" > "$OUT_DIR/${base}_consensus.fasta"
done

echo "✅ Tous les fichiers consensus ont été enregistrés dans : $OUT_DIR"
