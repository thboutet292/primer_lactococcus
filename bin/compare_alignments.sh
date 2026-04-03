#!/usr/bin/env bash
# =============================================================
# compare_alignments.sh
# ÉTAPE 4 du pipeline : Comparaison visuelle brut vs trimmé
#
# Ce script génère un graphique PNG à deux panneaux superposés
# qui permet de visualiser l'impact du nettoyage trimAl sur
# l'alignement multiple des séquences 16S de Lactococcus.
#
# Graphique 1 : Alignement MAFFT brut
#   - Ligne rouge  : % de gaps par colonne (bruit)
#   - Ligne bleue  : % de conservation par colonne (signal)
#   → Montre les "pics" de gaps qui perturbent l'analyse
#
# Graphique 2 : Alignement trimmé (post-trimAl)
#   - Même métriques mais sur l'alignement nettoyé
#   - Fond bleu rempli : met en valeur le signal biologique conservé
#   → Prouve que le nettoyage améliore la lisibilité du signal
#
# DÉPENDANCES :
#   pip install matplotlib
#
# SORTIE :
#   alignment/alignment_comparison.png 
# =============================================================

set -euo pipefail

log() { echo "[$(date +%H:%M:%S)] $*"; }

# Chemins des fichiers d'entrée
RAW_ALN="alignment/lactococcus_aln.fa"                # Alignement MAFFT brut
TRIMMED_ALN="alignment/lactococcus_aln_trimmed.fa"    # Alignement nettoyé
OUTPUT_PLOT="alignment/alignment_comparison.png"      # Fichier de sortie

log "Génération des graphiques de conservation en cours..."

python3 << PYEOF
import sys
import matplotlib.pyplot as plt
from collections import Counter

def parse_fasta_and_score(filepath):
    """
    Lit un fichier FASTA aligné et calcule pour chaque colonne :
        - gaps_pct  : pourcentage de gaps (mesure du bruit)
        - cons_pct  : pourcentage de conservation du nucléotide majoritaire
                      (calculé par rapport au nombre total de séquences,
                       les gaps étant donc pénalisés dans ce calcul)
    
    Paramètres :
        filepath : chemin vers le fichier FASTA aligné
    
    Retourne :
        gaps_pct (list), cons_pct (list), aln_len (int)
    """
    seqs = []
    curr_seq = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    # Nouvelle entrée : sauvegarde la séquence précédente
                    if curr_seq: seqs.append("".join(curr_seq))
                    curr_seq = []
                else:
                    curr_seq.append(line.upper())
            # Traitement de la dernière séquence
            if curr_seq: seqs.append("".join(curr_seq))
    except FileNotFoundError:
        print(f"Erreur : fichier {filepath} introuvable.")
        sys.exit(1)

    n_seq = len(seqs)
    aln_len = len(seqs[0])
    
    gaps_pct = []
    cons_pct = []
    
    # ── Analyse colonne par colonne ──
    for i in range(aln_len):
        # Extrait tous les caractères à la position i (une par séquence)
        col = [s[i] for s in seqs]
        
        # ── Calcul du taux de gaps ──
        # "-" = gap d'alignement, "N" = base indéterminée (traitée comme gap)
        n_gaps = sum(1 for b in col if b in ('-', 'N'))
        gaps_pct.append((n_gaps / n_seq) * 100)
        
        # ── Calcul du taux de conservation ──
        # On ne compte que les vraies bases (hors gaps/N)
        bases = [b for b in col if b not in ('-', 'N')]
        if bases:
            # Trouve la base la plus fréquente et son comptage
            top_base_count = Counter(bases).most_common(1)[0][1]
            # Conservation relative au total des séquences (gaps compris)
            # → pénalise les colonnes avec beaucoup de gaps
            cons_pct.append((top_base_count / n_seq) * 100)
        else:
            cons_pct.append(0)  # Colonne 100% gaps → conservation nulle
            
    return gaps_pct, cons_pct, aln_len

# ── Extraction des profils de conservation ──
gaps_raw, cons_raw, len_raw = parse_fasta_and_score("$RAW_ALN")
gaps_trim, cons_trim, len_trim = parse_fasta_and_score("$TRIMMED_ALN")

# ── Création de la figure à deux panneaux ──
# sharey=True : même axe Y pour les deux panneaux → comparaison directe
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharey=True)
fig.suptitle("Impact du nettoyage (trimAl) sur l'alignement multiple des 16S de Lactococcus", 
             fontsize=16, fontweight='bold')

# ── Panneau du HAUT : Alignement brut ──
# Montre la situation initiale avec les pics de gaps perturbateurs
ax1.set_title(f"Avant nettoyage (MAFFT Brut) - {len_raw} positions", fontsize=12)
ax1.plot(gaps_raw, color='red', alpha=0.7, linewidth=1.5, label='% Gaps (Bruit)')
ax1.plot(cons_raw, color='royalblue', alpha=0.9, linewidth=1.5, label='% Conservation (Signal)')
ax1.set_ylabel("Pourcentage (%)")
ax1.legend(loc="upper right")
ax1.grid(True, linestyle='--', alpha=0.5)

# ── Panneau du BAS : Alignement nettoyé ──
# fill_between remplit la zone sous la courbe de conservation
# en bleu semi-transparent → met visuellement en évidence les régions
# conservées qui seront utilisées pour le design d'amorces
ax2.set_title(f"Après nettoyage (trimAl) - {len_trim} positions", fontsize=12)
ax2.fill_between(range(len_trim), cons_trim, color='royalblue', alpha=0.3)
ax2.plot(cons_trim, color='royalblue', alpha=0.9, linewidth=1.5, label='% Conservation (Signal)')
ax2.plot(gaps_trim, color='red', alpha=0.7, linewidth=1.5, label='% Gaps (Bruit)')
ax2.set_ylabel("Pourcentage (%)")
ax2.set_xlabel("Position dans l'alignement (bp)")
ax2.legend(loc="upper right")
ax2.grid(True, linestyle='--', alpha=0.5)

# ── Sauvegarde du graphique ──
# bbox_inches='tight' : évite que les labels soient coupés en bord de figure
# dpi=300 : résolution haute qualité pour publication ou rapport
plt.tight_layout()
plt.savefig("$OUTPUT_PLOT", dpi=300, bbox_inches='tight')
print(f"Graphique généré avec succès : $OUTPUT_PLOT")
PYEOF

log "Terminé."