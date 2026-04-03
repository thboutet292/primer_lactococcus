#!/usr/bin/env bash
# =============================================================
# lactococcus_align.sh
# ÉTAPE 2 du pipeline : Alignement multiple + score de conservation
#
# Ce script aligne toutes les séquences 16S de Lactococcus avec
# MAFFT, puis calcule pour chaque colonne de l'alignement un score
# de conservation (fraction de séquences partageant la même base).
#
# PRÉREQUIS :
#   - build_databases.sh doit avoir été exécuté
#   - databases/Lactococcus_nr.fa doit exister
#
# DÉPENDANCES :
#   - mafft   : alignement multiple progressif
#   - python3 : calcul des scores de conservation
#
# SORTIES (dans alignment/) :
#   lactococcus_aln.fa              ← alignement FASTA → étapes suivantes
#   lactococcus_aln.clw             ← alignement ClustalW → Jalview/Seaview
#   mafft.log                       ← log complet de MAFFT
# =============================================================

set -euo pipefail

# ─────────────────────────────────────────────
# CHEMINS (exportés car utilisés dans les blocs Python)
# ─────────────────────────────────────────────
export DBDIR="$HOME/primer_lactococcus/databases"   # Dossier source contenant l BDD de Lactococcus 
export OUTDIR="$HOME/primer_lactococcus/alignment"  # Dossier de sortie des fichiers d'alignements 

# ─────────────────────────────────────────────
# PARAMÈTRES BIOLOGIQUES
# ─────────────────────────────────────────────
THREADS=4       # Nombre de coeurs alloués à MAFFT
MIN_WINDOW=20   # En dessous de 20 pb, une amorce serait trop courte pour être spécifique
MIN_CONS=0.90   # Fraction minimale de séquences devant partager la même base en une colonne 
MAX_GAP=0.05    # 0.05 = au plus 5% des séquences peuvent avoir un gap en cette position
  
# ─────────────────────────────────────────────
# FONCTIONS UTILITAIRES
# ─────────────────────────────────────────────
log()        { echo "[$(date +%H:%M:%S)] $*"; }     # Date
die()        { echo "ERREUR : $*" >&2; exit 1; }    # Si erreur 

# Compte le nombre de séquences dans un fichier FASTA
count_seqs() { grep -c "^>" "$1" 2>/dev/null || echo 0; }

# ─────────────────────────────────────────────
# VÉRIFICATION DES DÉPENDANCES
# ─────────────────────────────────────────────
check_deps() {
    log "Verification des dependances..."
    for cmd in mafft python3; do
        command -v "$cmd" >/dev/null 2>&1 \
            || die "$cmd introuvable — verifie que l environnement conda est active"
    done

    # Vérifie que le fichier d'entrée principal est présent
    [[ -f "$DBDIR/Lactococcus_nr.fa" ]] \
        || die "Lactococcus_nr.fa introuvable — lance d abord 01_build_databases.sh"
    log "  -> OK"
}

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
check_deps
mkdir -p "$OUTDIR"
log "Dossier de sortie : $OUTDIR"
echo ""

# =============================================================
# Alignement multiple avec MAFFT
# =============================================================
# --auto    : MAFFT choisit automatiquement la stratégie optimale
#             selon la taille et la nature des séquences
# --reorder : réorganise les séquences par similarité dans la sortie
#
# On génère DEUX formats de sortie :
#   - ClustalW (.clw) : pour la visualisation dans Jalview/Seaview
#   - FASTA    (.fa)  : pour le traitement Python en aval

log "=== Alignement multiple MAFFT ==="
log "  $(count_seqs \"$DBDIR/Lactococcus_nr.fa\") sequences a aligner sur $THREADS threads..."

# Format ClustalW
mafft \
    --auto \
    --thread "$THREADS" \
    --reorder \
    --clustalout \
    "$DBDIR/Lactococcus_nr.fa" \
    > "$OUTDIR/lactococcus_aln.clw" \
    2> "$OUTDIR/mafft.log"

# Format FASTA
mafft \
    --auto \
    --thread "$THREADS" \
    --reorder \
    "$DBDIR/Lactococcus_nr.fa" \
    > "$OUTDIR/lactococcus_aln.fa" \
    2>> "$OUTDIR/mafft.log"

log "  -> lactococcus_aln.clw  (ClustalW, visualisable dans Jalview/Seaview)"
log "  -> lactococcus_aln.fa   (FASTA, utilise pour le scoring)"
log "  -> mafft.log            (log MAFFT)"
echo ""

# ─────────────────────────────────────────────
# RÉSUMÉ FINAL
# ─────────────────────────────────────────────
log "================================================================"
log "Fichiers generes dans $OUTDIR :"
echo ""
printf "  %-38s %s\n" "lactococcus_aln.clw"          "ClustalW         -> Jalview / Seaview"
printf "  %-38s %s\n" "lactococcus_aln.fa"            "FASTA aligne     -> etapes suivantes"
printf "  %-38s %s\n" "mafft.log"                     "Log MAFFT"
echo ""
log "================================================================"