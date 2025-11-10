#!/bin/bash

CONS_DIR="/home/abo/stage/dengue_virus/Illumina0/DENV/ALLDENV/consensus_results"
cd "$CONS_DIR" || { echo "❌ Impossible d'accéder au dossier."; exit 1; }

rm -f all_consensus.fasta aligned_consensus.fasta

echo "Concaténation des fichiers consensus..."
cat *_consensus.fasta > all_consensus.fasta
if [ ! -s all_consensus.fasta ]; then
  echo "❌ Le fichier all_consensus.fasta est vide."
  exit 1
fi

echo "Alignement MAFFT (1 thread, 100 itérations max)..."
mafft --retree 1 --maxiterate 100 --thread 1 all_consensus.fasta > aligned_consensus.fasta
if [ $? -ne 0 ] || [ ! -s aligned_consensus.fasta ]; then
  echo "❌ L'alignement a échoué."
  exit 1
fi

echo "✅ Alignement réussi : $CONS_DIR/aligned_consensus.fasta"

