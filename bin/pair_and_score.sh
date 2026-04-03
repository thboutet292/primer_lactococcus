#!/usr/bin/env bash
# =============================================================
# pair_and_score.sh
# ÉTAPE 7 du pipeline : Appariement asymétrique + scoring thermodynamique
#
# Ce script constitue la DERNIÈRE étape analytique du pipeline.
# Il assemble des paires d'amorces "asymétriques" :
#   - Forward UNIVERSELLE : générée de novo depuis les fenêtres conservées
#     (bonne couverture de tous les Lactococcus, pas forcément spécifique)
#   - Reverse ULTRA-SPÉCIFIQUE : issue du criblage check_specificity.sh
#     (garantit la spécificité du genre, le "verrou" est côté Reverse)
#
# STRATÉGIE ASYMÉTRIQUE :
#   Cette approche est un compromis intelligent : on accepte une Forward
#   semi-universelle (qui amplifiera de nombreux Lactobacillales) mais on
#   compte sur la Reverse hyper-spécifique pour que le produit PCR COMPLET
#   soit exclusif à Lactococcus. La Reverse "verrouille" la spécificité.
#
# CRITÈRES DE SÉLECTION DES PAIRES :
#   1. La Forward doit être EN AMONT (position < ) de la Reverse
#   2. La taille d'amplicon doit être : MIN_AMP <= taille <= MAX_AMP (pb)
#   3. Le ΔTm (différence de Tm entre Forward et Reverse) doit être <= 5°C
#      → Si les deux amorces ont des Tm trop différentes, l'une sera
#        sur-utilisée et l'autre sous-utilisée → amplification inefficace
#
# CALCULS THERMODYNAMIQUES UTILISÉS :
#   - Tm (Melting Temperature) : formule de Wallace modifiée
#       Tm = 64.9 + 41 × (nGC - 16.4) / longueur
#       (valide pour les amorces courtes, 14 nt minimum)
#   - %GC : fraction de G et C dans la séquence
#   - has_run : détecte les répétitions de 4 bases identiques
#               (ex: AAAA, CCCC) qui déstabilisent l'hybridation
#
# DÉPENDANCES :
#   - python3 (modules standard uniquement : csv, os, collections)
#
# SORTIES :
#   specificity/asymmetric_pairs.tsv        ← tableau de toutes les paires
#   specificity/asymmetric_pairs_report.txt ← top 20 paires lisibles
# =============================================================

set -euo pipefail

# Fichiers d'entrée
WINDOWS_FASTA="alignment/conserved_windows_trimmed.fa"            # Régions conservées
VALIDATED_REV="specificity/validated_ultra_specific_primers.tsv"  # Reverse validées
OUTDIR="specificity"                                              # Dossier de sortie

# Taille minimale et maximale de l'amplicon à produire (pb)
MIN_AMP=150
MAX_AMP=800

log() { echo "[$(date +%H:%M:%S)] $*"; }

log "=== Recherche des paires asymétriques idéales ==="

python3 << PYEOF
import csv, os
from collections import defaultdict

win_file = "$WINDOWS_FASTA"
rev_file = "$VALIDATED_REV"
out_tsv  = "$OUTDIR/asymmetric_pairs.tsv"
out_txt  = "$OUTDIR/asymmetric_pairs_report.txt"

# ─────────────────────────────────────────────
# FONCTIONS PHYSICO-CHIMIQUES
# ─────────────────────────────────────────────

def tm(seq):
    """
    Calcule la température de fusion (Tm) d'une amorce courte.
    
    Utilise la formule de Wallace modifiée, valide pour les
    oligonucléotides de 14 à ~30 nt.
    
        Tm = 64.9 + 41 × (nGC - 16.4) / longueur
    
    Pour les amorces < 14 nt, utilise la formule basique :
        Tm = 2 × nAT + 4 × nGC
    (car les formules basées sur la longueur sont moins fiables)
    
    Paramètres :
        seq : séquence nucléotidique (str)
    
    Retourne :
        float : Tm en °C
    """
    seq = seq.upper()
    gc  = sum(1 for b in seq if b in "GC")
    at  = sum(1 for b in seq if b in "AT")
    return 64.9 + 41 * (gc - 16.4) / len(seq) if len(seq) >= 14 else 2 * at + 4 * gc

def gc_pct(seq):
    """
    Calcule le pourcentage de G+C d'une amorce.
    Idéal pour des amorces PCR : entre 40% et 65%.
    
    Paramètres :
        seq : séquence nucléotidique (str)
    
    Retourne :
        float : pourcentage GC (0 à 100)
    """
    return sum(1 for b in seq.upper() if b in "GC") / len(seq) * 100

def has_run(seq, n=4):
    """
    Détecte les "runs" : répétitions de n bases identiques consécutives.
    
    Exemple : "AAAA", "CCCC", "TTTTG" → problématiques car :
    - Provoquent des structures secondaires intramoléculaires
    - Réduisent l'efficacité d'hybridation de l'amorce
    - Peuvent causer des glissements de l'ADN polymérase
    
    Paramètres :
        seq : séquence nucléotidique (str)
        n   : longueur minimale du run à détecter (défaut : 4)
    
    Retourne :
        bool : True si un run est détecté, False sinon
    """
    for base in "ACGT":
        if base * n in seq.upper(): return True
    return False

def reverse_complement(seq):
    """
    Calcule le complément inverse d'une séquence nucléotidique.
    
    Nécessaire pour les amorces Reverse : la séquence stockée dans
    les fichiers FASTA est en orientation 5'→3' sur le brin sens.
    Pour une amorce Reverse, on doit l'inverser et la complémenter
    pour obtenir sa séquence telle qu'on la synthétise (5'→3' antisens).
    
    Paramètres :
        seq : séquence ADN (str), peut contenir N et -
    
    Retourne :
        str : complément inverse en majuscules
    """
    comp = {"A":"T","T":"A","G":"C","C":"G", "N":"N","-":"-"}
    return "".join(comp.get(b,"N") for b in reversed(seq.upper()))

# ─────────────────────────────────────────────
# 1. CHARGEMENT DES FENÊTRES ET GÉNÉRATION DES AMORCES FORWARD
# ─────────────────────────────────────────────
# Pour chaque fenêtre conservée, on génère des sous-séquences de
# 18 à 22 nt avec une fenêtre glissante. Seules celles qui passent
# les filtres thermodynamiques sont gardées comme candidates Forward.
fwd_candidates = []
with open(win_file, "r") as f:
    win_id    = ""
    win_start = 0
    seq_acc   = []
    
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # ── Traitement de la fenêtre précédente ──
            if win_id and "".join(seq_acc):
                seq = "".join(seq_acc)
                # Sous-séquences de 18 à 22 nt (toutes les tailles d'amorces)
                for size in range(18, 23):
                    for i in range(len(seq) - size + 1):
                        subseq = seq[i:i+size]
                        t, g = tm(subseq), gc_pct(subseq)
                        # Filtres thermodynamiques :
                        #   45°C <= Tm <= 65°C  : plage d'hybridation PCR standard
                        #   40%  <= GC <= 65%   : stabilité correcte sans sur-stabilité
                        #   has_run = False      : pas de répétitions problématiques
                        if 45 <= t <= 65 and 40 <= g <= 65 and not has_run(subseq):
                            fwd_candidates.append({
                                "id":           f"{win_id}_Fwd_pos{i}",
                                "seq":          subseq,
                                "win":          win_id,
                                "start_global": win_start + i,  # Position dans l'alignement
                                "len":          size,
                                "tm":           round(t, 1),
                                "gc":           round(g, 1)
                            })
            
            # ── Parsing de l'en-tête FASTA ──
            # Format attendu : >W001_pos123-456_len333_cons0.97
            parts      = line.split("_")
            win_id     = parts[0][1:]  # "W001" (sans le ">")
            # Extrait la position de début depuis "pos123-456"
            win_start  = int(parts[1].replace("pos", "").split("-")[0])
            seq_acc    = []
        else:
            seq_acc.append(line.upper())
            
    # ── Traitement de la dernière fenêtre ──
    if win_id and "".join(seq_acc):
        seq = "".join(seq_acc)
        for size in range(18, 23):
            for i in range(len(seq) - size + 1):
                subseq = seq[i:i+size]
                t, g = tm(subseq), gc_pct(subseq)
                if 45 <= t <= 65 and 40 <= g <= 65 and not has_run(subseq):
                    fwd_candidates.append({
                        "id": f"{win_id}_Fwd_pos{i}",
                        "seq": subseq, "win": win_id,
                        "start_global": win_start + i, "len": size,
                        "tm": round(t, 1), "gc": round(g, 1)
                    })

print(f"  {len(fwd_candidates)} amorces Forward potentielles générées depuis l'alignement.")

# ─────────────────────────────────────────────
# 2. CHARGEMENT DES AMORCES REVERSE VALIDÉES
# ─────────────────────────────────────────────
# On charge les amorces Reverse produites par check_specificity.sh
# et on reconstitue leur position globale dans l'alignement
# en relisant les en-têtes du fichier de fenêtres
rev_candidates = []
with open(rev_file, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        pid = row["Primer_ID"]
        seq = row["Sequence_5_to_3"]
        win = pid.split("_")[0]  # Fenêtre d'origine (ex: "W003")
        # Position de la Reverse dans sa fenêtre (ex: "pos5" → 5)
        pos_in_win = int(pid.split("_")[-1].replace("pos", ""))
        
        # ── Recherche de la position globale de la fenêtre ──
        # On relit le fichier de fenêtres pour retrouver le start_global
        start_global = 0
        with open(win_file, "r") as wf:
            for line in wf:
                if line.startswith(f">{win}"):
                    start_global = int(line.strip().split("_")[1].replace("pos", "").split("-")[0])
                    break
                    
        rev_candidates.append({
            "id":           pid,
            "seq":          seq,
            "win":          win,
            "start_global": start_global + pos_in_win,  # Position absolue dans l'alignement
            "len":          len(seq),
            "tm":           round(tm(seq), 1),
            "gc":           round(gc_pct(seq), 1)
        })

print(f"  {len(rev_candidates)} amorces Reverse Ultra-Spécifiques chargées.")

# ─────────────────────────────────────────────
# 3. APPARIEMENT FORWARD × REVERSE
# ─────────────────────────────────────────────
# Pour chaque couple (Forward, Reverse), on vérifie :
#   1. La Forward est bien EN AMONT (position < ) de la Reverse
#      (sinon on ne peut pas amplifier entre les deux)
#   2. La taille d'amplicon estimée est dans la plage acceptable
#   3. Le ΔTm est <= 5°C (compatibilité de PCR)
#
# CALCUL DE LA TAILLE D'AMPLICON :
#   amplicon = (fin de la Reverse) - (début de la Forward)
#            = (start_global_rev + len_rev) - start_global_fwd
#   Cette estimation est basée sur les positions dans l'alignement.
pairs = []
for fwd in fwd_candidates:
    for rev in rev_candidates:
        # ── Test 1 : orientation (Forward avant Reverse) ──
        if fwd["start_global"] >= rev["start_global"]:
            continue  # La Forward est en aval de la Reverse → invalide
            
        # ── Test 2 : taille d'amplicon ──
        amp_size = (rev["start_global"] + rev["len"]) - fwd["start_global"]
        
        if not ($MIN_AMP <= amp_size <= $MAX_AMP):
            continue  # Amplicon trop court ou trop long
            
        # ── Test 3 : compatibilité thermodynamique ──
        delta_tm = abs(fwd["tm"] - rev["tm"])
        if delta_tm > 5:
            continue  # Trop grande différence de Tm → PCR inefficace
            
        # La paire passe tous les critères : on l'enregistre
        pairs.append({
            "fwd_id":   fwd["id"],
            "fwd_seq":  fwd["seq"],
            "fwd_tm":   fwd["tm"],
            "fwd_gc":   fwd["gc"],
            "rev_id":   rev["id"],
            # La Reverse est stockée en complément inverse (orientation synthèse)
            "rev_seq":  reverse_complement(rev["seq"]),
            "rev_tm":   rev["tm"],
            "rev_gc":   rev["gc"],
            "delta_tm": round(delta_tm, 1),
            "amp_size": amp_size,
            "win_fwd":  fwd["win"],
            "win_rev":  rev["win"]
        })

print(f"  {len(pairs)} paires asymétriques viables trouvées !")

# ─────────────────────────────────────────────
# 4. EXPORT DES RÉSULTATS
# ─────────────────────────────────────────────
# Les paires sont triées par ΔTm croissant (meilleure compatibilité en tête)
# puis par taille d'amplicon décroissante (amplicons plus longs préférés
# pour une meilleure résolution électrophorétique)

# ── Export TSV (tableau complet de toutes les paires) ──
fields = ["fwd_id","fwd_seq","fwd_tm","fwd_gc","rev_id","rev_seq","rev_tm","rev_gc","delta_tm","win_fwd","win_rev","amp_size"]
with open(out_tsv, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for p in sorted(pairs, key=lambda x: (x["delta_tm"], -x["amp_size"])):
        w.writerow(p)

# ── Export rapport lisible (top 20 paires) ──
with open(out_txt, "w") as out:
    out.write("=" * 80 + "\n")
    out.write(f"  TOP PAIRES ASYMÉTRIQUES (Forward Universelle x Reverse Spécifique)\n")
    out.write("=" * 80 + "\n\n")
    for i, p in enumerate(sorted(pairs, key=lambda x: (x["delta_tm"], -x["amp_size"]))[:20], 1):
        out.write(f"Paire #{i:02d}  [{p['win_fwd']} x {p['win_rev']}]  |  Amplicon estimé : ~{p['amp_size']} pb\n")
        out.write(f"  Forward : {p['fwd_seq']:<25} Tm={p['fwd_tm']}°C  GC={p['fwd_gc']}%\n")
        out.write(f"  Reverse : {p['rev_seq']:<25} Tm={p['rev_tm']}°C  GC={p['rev_gc']}%\n")
        out.write(f"  DeltaTm : {p['delta_tm']}°C\n\n")

PYEOF

log "================================================================"
cat $OUTDIR/asymmetric_pairs_report.txt
log "================================================================"