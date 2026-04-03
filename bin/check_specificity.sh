#!/usr/bin/env bash
# =============================================================
# check_specificity.sh
# ÉTAPE 6 du pipeline : Criblage de spécificité in silico
#
# Ce script teste la spécificité de toutes les amorces candidates
# contre une base de séquences "outgroups" (genres proches de
# Lactococcus). Une amorce est validée UNIQUEMENT si elle ne peut
# pas s'hybrider sur les séquences des genres non-cibles.
#
# STRATÉGIE DE VALIDATION ("Ultra-Specific 3' lock") :
#   1. Toutes les sous-séquences de 18 à 22 nt des fenêtres
#      conservées sont testées comme amorces candidates.
#   2. Pour chaque candidate, on la compare à chaque sous-séquence
#      de même taille dans tous les outgroups.
#   3. Une candidate est REJETÉE si un outgroup présente :
#         - <= MAX_MISMATCHES_ALLOWED mismatches sur toute sa longueur
#      ET - 0 mismatch sur les CRITICAL_3P_LEN derniers nucléotides (3')
#      
#      POURQUOI LE CRITÈRE 3' EST-IL CRITIQUE ?
#      L'extrémité 3' d'une amorce est le point d'amorçage de
#      l'ADN polymérase. Si les dernières bases (5 nt) s'hybrident
#      parfaitement sur un outgroup, la polymérase peut s'étendre
#      même si quelques mismatches existent en 5'. Une amorce qui
#      ferait cela serait non-spécifique.
#      Le "3' lock" garantit qu'au moins un mismatch existe en 3'
#      sur TOUS les outgroups → l'hybridation non-spécifique est bloquée.
#
# DÉPENDANCES :
#   - python3 (modules standard uniquement : os, sys, collections)
#
# SORTIE :
#   specificity/validated_ultra_specific_primers.tsv
#   Colonnes : Primer_ID | Sequence_5_to_3 | Length | Status | Comment
# =============================================================

set -euo pipefail

echo "[$(date +%H:%M:%S)] Démarrage du cribleur absolu (Version Dynamique & Transparente)..."

python3 << 'PYEOF'
import os
import sys
from collections import Counter

# ─────────────────────────────────────────────
# PARAMÈTRES BIOLOGIQUES ET SEUILS
# ─────────────────────────────────────────────

# Longueur minimale d'une amorce testée (nt)
# En dessous de 18 nt → risque de manque de spécificité
MIN_PRIMER_LEN = 18

# Longueur maximale d'une amorce testée (nt)
# Au-delà de 22 nt → coût de synthèse plus élevé, avantage limité
MAX_PRIMER_LEN = 22

# Nombre maximum de mismatches tolérés sur la longueur totale de l'amorce
# 3 mismatches = la candidate peut avoir jusqu'à 3 différences avec un outgroup
# sans déclencher de rejet (seul le critère 3' est éliminatoire)
MAX_MISMATCHES_ALLOWED = 3

# Longueur de l'extrémité 3' soumise au test de spécificité strict
# 5 nt = les 5 dernières bases doivent TOUJOURS avoir >= 1 mismatch sur les outgroups
CRITICAL_3P_LEN = 5

# Fichiers d'entrée
WINDOWS_FASTA  = "alignment/conserved_windows_trimmed.fa"  # Régions candidates
OUTGROUPS_FASTA = "databases/outgroups_all_nr.fa"           # Séquences à éviter

# Fichier de sortie
OUTPUT_TSV = "specificity/validated_ultra_specific_primers.tsv"

def read_fasta(filepath):
    """
    Lit un fichier FASTA et retourne un dictionnaire {id_sequence: séquence}.
    
    Paramètres :
        filepath : chemin vers le fichier FASTA
    
    Retourne :
        dict {str: str} — clés = identifiants, valeurs = séquences en majuscules
    """
    seqs = {}
    curr_id, curr_seq = "", []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if curr_id: seqs[curr_id] = "".join(curr_seq).upper()
                curr_id = line[1:].split()[0]  # Prend uniquement le premier mot
                curr_seq = []
            else:
                curr_seq.append(line.upper())
        if curr_id: seqs[curr_id] = "".join(curr_seq).upper()
    return seqs

print("  Chargement des séquences en mémoire...")
try:
    windows   = read_fasta(WINDOWS_FASTA)   # Fenêtres conservées (peu nombreuses)
    outgroups = read_fasta(OUTGROUPS_FASTA) # Outgroups (plusieurs milliers de séquences)
except FileNotFoundError as e:
    print(f"Erreur : {e}")
    sys.exit(1)

# Crée le dossier de sortie si inexistant
os.makedirs(os.path.dirname(OUTPUT_TSV), exist_ok=True)

valid_primers_count = 0
rejection_reasons = Counter()  # Compteur : quel outgroup rejette le plus d'amorces

with open(OUTPUT_TSV, "w") as out_file:
    out_file.write("Primer_ID\tSequence_5_to_3\tLength\tStatus\tComment\n")
    
    # ── BOUCLE PRINCIPALE : fenêtres conservées ──
    for w_name, w_seq in windows.items():
        # Préfixe court pour l'identifiant (ex: "W001" depuis "W001_pos123-456")
        short_w = w_name.split("_")[0]
        
        # ── TEST DE TOUTES LES TAILLES D'AMORCES DANS CHAQUE FENÊTRE ──
        # On génère des sous-séquences ("amorces candidates") de 18 à 22 nt
        # en fenêtre glissante sur la séquence consensus de la fenêtre
        for p_len in range(MIN_PRIMER_LEN, MAX_PRIMER_LEN + 1):
            if len(w_seq) < p_len: continue  # Fenêtre trop courte pour cette taille
            
            # ── Fenêtre glissante sur la séquence consensus ──
            for i in range(len(w_seq) - p_len + 1):
                primer_cand = w_seq[i : i + p_len]  # Amorce candidate
                primer_3p   = primer_cand[-CRITICAL_3P_LEN:]  # Extrémité 3' critique
                
                is_safe    = True    # Hypothèse : l'amorce est spécifique
                killer_taxon = None  # Quel outgroup la rejette (pour le rapport)
                
                # ── CRIBLAGE CONTRE TOUS LES OUTGROUPS ──
                for out_id, out_seq in outgroups.items():
                    if not is_safe: break  # Court-circuit : déjà rejetée
                    
                    # Fenêtre glissante sur la séquence outgroup
                    for j in range(len(out_seq) - p_len + 1):
                        target = out_seq[j : j + p_len]
                        
                        # ── TEST 1 : mismatches sur la longueur totale ──
                        # Compare base par base l'amorce et la cible outgroup
                        mismatches = sum(1 for a, b in zip(primer_cand, target) if a != b)
                        
                        # Si trop similaire globalement, on inspecte le 3'
                        if mismatches <= MAX_MISMATCHES_ALLOWED:
                            target_3p   = target[-CRITICAL_3P_LEN:]
                            mismatches_3p = sum(1 for a, b in zip(primer_3p, target_3p) if a != b)
                            
                            # ── TEST 2 : mismatch obligatoire en 3' ──
                            # Si l'extrémité 3' s'hybride parfaitement (0 mismatch)
                            # sur un outgroup → l'ADN polymérase pourrait s'étendre
                            # → l'amorce est NON-SPÉCIFIQUE → REJET
                            if mismatches_3p == 0:
                                is_safe      = False
                                # Extrait le genre outgroup (1er segment avant "|")
                                # Ex: "Streptococcus|RefSeq|NR_044534" → "Streptococcus"
                                killer_taxon = out_id.split('|')[0]
                                break  # Inutile de continuer : l'amorce est rejetée
                
                # ── BILAN POUR CETTE AMORCE ──
                if is_safe:
                    # L'amorce passe tous les tests : on l'enregistre
                    primer_id = f"{short_w}_len{p_len}_pos{i}"
                    out_file.write(f"{primer_id}\t{primer_cand}\t{p_len}\tVALIDATED\tUltra-specific 3' lock\n")
                    valid_primers_count += 1
                else:
                    # L'amorce est rejetée : on note le genre responsable
                    rejection_reasons[killer_taxon] += 1

# ── RAPPORT FINAL ──
print("\n" + "="*60)
print("             RAPPORT DE SPÉCIFICITÉ IN SILICO")
print("="*60)
print(f"Amorces testées (tailles 18 à 22 nt) : {sum(rejection_reasons.values()) + valid_primers_count}")
print(f"Amorces VALIDES (Ultra-Spécifiques)  : {valid_primers_count}")
print("\nPourquoi les amorces ont-elles échoué ? (Top des faux-positifs) :")

# Affiche les 5 genres outgroup qui ont rejeté le plus d'amorces
for taxon, count in rejection_reasons.most_common(5):
    print(f"Bloqué par {taxon:<20} : {count} amorces rejetées")

print("="*60)
if valid_primers_count == 0:
    # Aucune amorce n'a passé le criblage → diagnostic biologique
    print("CONCLUSION BIOLOGIQUE :")
    print("Le gène 16S ne possède pas de variabilité suffisante en 3'")
    print("pour différencier Lactococcus des taxons listés ci-dessus.")
print(f"-> Résultats sauvegardés dans : {OUTPUT_TSV}\n")
PYEOF