#!/usr/bin/env bash
# =============================================================
# extract_window.sh
# ÉTAPE 5 du pipeline : Extraction des régions consensus conservées
#
# Ce script identifie dans l'alignement trimmé les "fenêtres
# conservées" : des régions contiguës où la même base est
# majoritairement partagée par toutes les séquences de Lactococcus.
# La séquence consensus de chaque fenêtre est ensuite exportée en
# FASTA et servira de point de départ pour générer les amorces
# candidates dans les étapes suivantes.
#
# DIFFÉRENCE AVEC lactococcus_trim.sh (étape 3) :
#   lactococcus_trim.sh identifie aussi des fenêtres, mais c'est
#   une étape intermédiaire de validation. extract_window.sh est
#   l'étape officielle de production qui génère le fichier
#   conserved_windows_trimmed.fa utilisé par check_specificity.sh.
#   Les paramètres (MIN_CONS, MAX_GAP) peuvent être ajustés ici
#   indépendamment pour explorer différents jeux de critères.
#
# ALGORITHME EN DEUX TEMPS :
#   1. Pour chaque colonne de l'alignement :
#      → Calcule le taux de gaps et le taux de conservation
#      → Marque la colonne "bonne" si elle passe les deux seuils
#   2. Fusionne les colonnes "bonnes" adjacentes en blocs contigus
#      → Un bloc de longueur >= MIN_WINDOW devient une fenêtre
#      → Extrait la base consensus (majoritaire) de chaque colonne
#         pour construire la séquence consensus de la fenêtre
#
# DÉPENDANCES :
#   - python3 (modules standard uniquement : os, sys, collections)
#
# SORTIES :
#   alignment/conserved_windows_trimmed.fa  ← entrée de check_specificity.sh
# =============================================================

set -euo pipefail

log() { echo "[$(date +%H:%M:%S)] $*"; }

# ─────────────────────────────────────────────
# CONFIGURATION ET PARAMÈTRES
# ─────────────────────────────────────────────

export INPUT_ALN="alignment/lactococcus_aln_trimmed.fa"         # Fichier d'entrée alignement trimmé
export OUTPUT_FASTA="alignment/conserved_windows_trimmed.fa"    # Fichier de sortie séquences consensus des fenêtres conservées

export MIN_WINDOW=18                                            # En dessous de 18 pb, une amorce serait trop courte pour être spécifique
export MIN_CONS=0.70                                            # 0.70 = au moins 70% des séquences doivent avoir la même base
export MAX_GAP=0.05                                             # 0.05 = au plus 5% des séquences peuvent avoir un gap à cette position

log "Recherche des fenêtres conservées sur : $INPUT_ALN"
log "Paramètres : MIN_WINDOW=$MIN_WINDOW, MIN_CONS=$MIN_CONS, MAX_GAP=$MAX_GAP"

# ─────────────────────────────────────────────
# SCRIPT PYTHON D'EXTRACTION
# ─────────────────────────────────────────────
python3 << 'PYEOF'
import os
import sys
from collections import Counter

def run_extraction():
    """
    Fonction principale : identifie et exporte les fenêtres conservées
    de l'alignement trimmé selon les critères de conservation et de gaps.
    """
    # ── Lecture des paramètres depuis les variables d'environnement bash ──
    input_file  = os.environ.get("INPUT_ALN")
    output_file = os.environ.get("OUTPUT_FASTA")
    min_w = int(os.environ.get("MIN_WINDOW"))     # Longueur minimale d'une fenêtre
    min_c = float(os.environ.get("MIN_CONS"))     # Conservation minimale
    max_g = float(os.environ.get("MAX_GAP"))      # Taux de gaps maximal

    # ── Lecture du fichier FASTA aligné ──
    seqs = []       # Liste de séquences (une chaîne par séquence)
    curr_seq = []   # Accumulateur de lignes pour la séquence courante

    try:
        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    # Nouvelle séquence : sauvegarde la précédente
                    if curr_seq:
                        seqs.append("".join(curr_seq))
                    curr_seq = []
                else:
                    curr_seq.append(line.upper())
            # Traitement de la dernière séquence
            if curr_seq:
                seqs.append("".join(curr_seq))
    except FileNotFoundError:
        print(f"Erreur : {input_file} introuvable.")
        sys.exit(1)

    if not seqs:
        print("Erreur : Alignement vide.")
        sys.exit(1)

    n_seq   = len(seqs)        # Nombre total de séquences dans l'alignement
    aln_len = len(seqs[0])    # Longueur de l'alignement (colonnes)

    # ── PREMIÈRE PASSE : identification des colonnes "bonnes" ──
    # Une colonne est "bonne" si :
    #   1. Son taux de gaps est <= MAX_GAP
    #   2. Sa conservation (base majoritaire / total séquences) est >= MIN_CONS
    good_positions = []
    for i in range(aln_len):
        # Extrait la colonne i : un caractère par séquence
        col = [s[i] for s in seqs]
        
        # ── Test 1 : taux de gaps ──
        n_gaps   = sum(1 for b in col if b in ('-', 'N'))
        gap_frac = n_gaps / n_seq
        
        if gap_frac <= max_g:
            # ── Test 2 : conservation de la base majoritaire ──
            # On exclut les gaps/N avant de chercher la base majoritaire
            bases_only = [b for b in col if b not in ('-', 'N')]
            if bases_only:
                top_base, top_count = Counter(bases_only).most_common(1)[0]
                # Conservation = nb occurrences de la base majoritaire / nb total séquences
                # (n_seq au dénominateur, pas len(bases_only) → les gaps pénalisent)
                if (top_count / n_seq) >= min_c:
                    good_positions.append(i)  # Cette colonne passe les deux tests

    # ── DEUXIÈME PASSE : fusion en fenêtres contiguës ──
    # On parcourt les positions "bonnes" triées et on fusionne les adjacentes
    windows = []
    if good_positions:
        start = prev = good_positions[0]
        for pos in good_positions[1:]:
            if pos == prev + 1:
                prev = pos  # Position consécutive → on étend la fenêtre courante
            else:
                # Discontinuité → on ferme la fenêtre si elle atteint la taille minimale
                if (prev - start + 1) >= min_w:
                    windows.append((start, prev))
                start = prev = pos  # On commence une nouvelle fenêtre
        # Traitement de la dernière fenêtre (non fermée dans la boucle)
        if (prev - start + 1) >= min_w:
            windows.append((start, prev))

    # ── Export FASTA des séquences consensus ──
    # Pour chaque fenêtre, on reconstruit la séquence consensus :
    # à chaque position, on prend la base la plus fréquente (hors gaps/N)
    with open(output_file, 'w') as out:
        for i, (start, end) in enumerate(windows, 1):
            length    = end - start + 1
            consensus = ""

            for p in range(start, end + 1):
                col_bases = [s[p] for s in seqs if s[p] not in ('-', 'N')]
                # Base consensus : la plus représentée à cette position
                consensus += Counter(col_bases).most_common(1)[0][0]
            
            # En-tête FASTA : encode les coordonnées et la longueur
            # Format : >W001_pos123-456_len333
            header = f">W{i:03d}_pos{start}-{end}_len{length}"
            out.write(f"{header}\n{consensus}\n")
            
    print(f"Extraction terminée : {len(windows)} fenêtres trouvées.")
    print(f"Fichier de sortie : {output_file}")

run_extraction()
PYEOF

log "Processus terminé."