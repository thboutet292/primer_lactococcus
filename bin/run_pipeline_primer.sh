#!/usr/bin/env bash
# =============================================================
# run_pipeline_primer.sh
# Exécution séquentielle du pipeline de design d'amorces 16S
# =============================================================

# Arrêter le script immédiatement si une commande échoue
set -euo pipefail

log() { echo -e "\n\033[1;32m[$(date +%H:%M:%S)] === $1 ===\033[0m"; }

log "DÉMARRAGE DU PIPELINE COMPLET - LACTOCOCCUS"

# Étape 1 : Construction et nettoyage des bases de données
log "ETAPE 1 : build_databases.sh"
bash bin/build_databases.sh

# Étape 2 : Alignement multiple
log "ETAPE 2 : lactococcus_align.sh"
bash bin/lactococcus_align.sh

# Étape 3 : Nettoyage de l'alignement
log "ETAPE 3 : lactococcus_trim.sh"
bash bin/lactococcus_trim.sh

# Étape 4 : Comparaison visuelle des alignements
log "ETAPE 4 : compare_alignments.sh"
bash bin/compare_alignments.sh

# Étape 5 : Extraction des régions conservées
log "ETAPE 5 : extract_window.sh"
bash bin/extract_window.sh

# Étape 6 : Cribleur de spécificité in silico
log "ETAPE 6 : check_specificity.sh"
bash bin/check_specificity.sh

# Étape 7 : Appariement asymétrique et scoring thermodynamique
log "ETAPE 7 : pair_and_score.sh (Recherche de paires asymétriques)"
bash bin/pair_and_score.sh

log "PIPELINE TERMINÉ AVEC SUCCÈS !"
echo "-> Les meilleures paires d'amorces se trouvent dans : specificity/asymmetric_pairs_report.txt"