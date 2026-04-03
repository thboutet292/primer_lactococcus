#!/usr/bin/env bash
# =============================================================
# lactococcus_trim.sh
# ÉTAPE 3 du pipeline : Nettoyage de l'alignement + re-scoring
#
# POURQUOI NETTOYER L'ALIGNEMENT ?
#   Un alignement MAFFT brut contient de nombreuses colonnes
#   artificielles : positions où la majorité des séquences ont un
#   gap ("-") dû à des insertions rares ou à des artefacts de
#   séquençage. Ces colonnes constituent du "bruit" qui :
#     - gonfle artificiellement la longueur de l'alignement
#     - dilue le signal de conservation des vraies régions conservées
#     - complique l'identification des fenêtres candidates
#
#   trimAl élimine ces colonnes problématiques selon un seuil de
#   gap configurable, puis on recalcule proprement les scores sur
#   l'alignement nettoyé.
#
# CE QUE FAIT CE SCRIPT (3 étapes) :
#   1. Trimming avec trimAl (suppression des colonnes à fort taux de gaps)
#   2. Validation du trimming (comparaison stats avant/après)
#   3. Re-calcul des scores de conservation colonne par colonne
#
# DÉPENDANCES :
#   - trimal  : nettoyage de l'alignement multiple
#   - python3 : re-scoring et ré-identification des fenêtres
#
# SORTIES (dans alignment/) :
#   lactococcus_aln_trimmed.fa      ← alignement nettoyé → étapes suivantes
#   lactococcus_aln_trimmed.clw     ← version ClustalW → Jalview/Seaview
# =============================================================

set -euo pipefail

# ─────────────────────────────────────────────
# CHEMINS
# ─────────────────────────────────────────────
# Dossier contenant l'alignement produit par lactococcus_align.sh
ALNDIR="$HOME/primer_lactococcus/alignment"

# ─────────────────────────────────────────────
# PARAMÈTRES TRIMAL
# ─────────────────────────────────────────────
GT=0.20          # on conserve les colonnes où au moins 80% des séquences n'ont pas de gaps 

# ─────────────────────────────────────────────
# PARAMÈTRES RE-SCORING (après trimming)
# ─────────────────────────────────────────────

MIN_WINDOW=20    # Longueur minimale d'une fenêtre conservée (pb)
MIN_CONS=0.90    # Consensus : 90% des séquences portent la même base
MAX_GAP=0.05     # Fraction maximale de gaps tolérée par colonne après trimming

# ─────────────────────────────────────────────
# FONCTIONS UTILITAIRES
# ─────────────────────────────────────────────
log() { echo "[$(date +%H:%M:%S)] $*"; }
die() { echo "ERREUR : $*" >&2; exit 1; }

# ─────────────────────────────────────────────
# VÉRIFICATION DES DÉPENDANCES
# ─────────────────────────────────────────────
check_deps() {
    log "Vérification des dépendances..."
    for cmd in trimal python3; do
        command -v "$cmd" >/dev/null 2>&1 \
            || die "$cmd introuvable — vérife que l'environnement conda est activé"
    done
    # Vérifie que le fichier d'entrée existe
    [[ -f "$ALNDIR/lactococcus_aln.fa" ]] \
        || die "lactococcus_aln.fa introuvable — lance d'abord lactococcus_align.sh"
    # Vérifie qu'il n'est pas vide (taille > 0)
    [[ -s "$ALNDIR/lactococcus_aln.fa" ]] \
        || die "lactococcus_aln.fa est vide — MAFFT a dû planter, relance lactococcus_align.sh"
    log "  → OK"
}

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
check_deps
echo ""

# =============================================================
# ÉTAPE 1/2 : Nettoyage avec trimAl
# =============================================================
# trimAl supprime les colonnes de l'alignement qui dépassent le
# seuil de gaps (-gt). On génère deux formats :
#   - FASTA  : utilisé pour le re-scoring Python ci-dessous
#   - ClustalW : pour la visualisation dans Jalview ou SeaView
log "=== ETAPE 1/2 : Trimming de l'alignement (trimal -gt $GT) ==="
log "  Fichier source : lactococcus_aln.fa"

# Format FASTA trimmé → utilisé pour le re-scoring Python
trimal \
    -in    "$ALNDIR/lactococcus_aln.fa" \
    -out   "$ALNDIR/lactococcus_aln_trimmed.fa" \
    -gt    "$GT" \
    -fasta

# Format ClustalW trimmé → visualisation Jalview / Seaview
trimal \
    -in    "$ALNDIR/lactococcus_aln.fa" \
    -out   "$ALNDIR/lactococcus_aln_trimmed.clw" \
    -gt    "$GT" \
    -clustal

log "  → lactococcus_aln_trimmed.fa   (FASTA, entrée re-scoring)"
log "  → lactococcus_aln_trimmed.clw  (ClustalW, visualisation Jalview/Seaview)"
echo ""

# =============================================================
# ÉTAPE 2/2 : Validation du trimming
# =============================================================
# Compare les statistiques de l'alignement avant et après trimming
# pour confirmer que les colonnes indésirables ont été supprimées
# sans perte du signal biologique principal.
log "=== ETAPE 2/2 : Validation du trimming ==="

# Note : le heredoc avec 'PYEOF' (apostrophes) désactive l'interpolation
# bash → les variables Python ($0, $i, etc.) ne sont pas interprétées par bash
python3 << 'PYEOF'
import os
from collections import Counter

def aln_stats(path, label):
    """
    Calcule et affiche les statistiques principales d'un fichier FASTA aligné.
    
    Paramètres :
        path  : chemin vers le fichier FASTA
        label : nom affiché dans la sortie (ex: "AVANT trimming")
    
    Retourne :
        la longueur de l'alignement (nombre de colonnes)
    """
    # ── Lecture du FASTA ──
    seqs, seq = [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if seq: seqs.append("".join(seq))
                seq = []
            else:
                seq.append(line.upper())
        if seq: seqs.append("".join(seq))

    n_seq   = len(seqs)
    aln_len = len(seqs[0]) if seqs else 0

    # ── Calcul du taux de gaps par colonne ──
    gap_fracs = []
    for i in range(aln_len):
        col = [s[i] for s in seqs]
        gap_fracs.append(sum(1 for b in col if b in ("-","N")) / n_seq)

    # Colonnes "propres" : au plus 5% de gaps
    cols_clean = sum(1 for g in gap_fracs if g <= 0.05)
    # Colonnes "très gappées" : plus de 50% de gaps → à supprimer
    cols_gappy = sum(1 for g in gap_fracs if g >  0.50)
    # Taux de gaps moyen sur tout l'alignement
    mean_gap   = sum(gap_fracs) / aln_len if aln_len else 0

    print(f"  [{label}]")
    print(f"    Séquences        : {n_seq}")
    print(f"    Longueur (col)   : {aln_len}")
    print(f"    Gap moyen/col    : {mean_gap:.3f}")
    print(f"    Colonnes gap>50% : {cols_gappy}")
    print(f"    Colonnes gap<=5% : {cols_clean}")
    return aln_len

alndir = os.path.expanduser("~/primer_lactococcus/alignment")
l_orig = aln_stats(os.path.join(alndir, "lactococcus_aln.fa"),         "AVANT trimming")
print()
l_trim = aln_stats(os.path.join(alndir, "lactococcus_aln_trimmed.fa"), "APRÈS trimming")
print()
removed = l_orig - l_trim
print(f"  → Colonnes supprimées : {removed} ({removed/l_orig*100:.1f}% de l'alignement original)")
PYEOF

echo ""



# ─────────────────────────────────────────────
# RÉSUMÉ FINAL
# ─────────────────────────────────────────────
log "================================================================"
log "Fichiers générés dans $ALNDIR :"
echo ""
printf "  %-42s %s\n" "lactococcus_aln_trimmed.fa"       "FASTA trimmé        → entrée check_specificity.sh"
printf "  %-42s %s\n" "lactococcus_aln_trimmed.clw"      "ClustalW trimmé     → Jalview / Seaview"
echo ""
log "================================================================"