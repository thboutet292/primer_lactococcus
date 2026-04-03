#!/usr/bin/env bash
# =============================================================
# build_databases.sh  (anciennement build_and_preprocess.sh)
# Pipeline COMPLET de préparation des bases de données 16S.
#   1. La BASE CIBLE (Lactococcus) :
#      Séquences 16S complètes (>= 1200 pb), dérepliquées, filtrées
#      des sous-séquences incluses et des chimères.
#      → Fichier final : databases/Lactococcus_nr.fa
#
#   2. Les BASES OUTGROUPS (6 genres proches) :
#      Séquences 16S (>= 900 pb), dérepliquées uniquement.
#      → Utilisées pour tester la spécificité des amorces candidates.
#      → Fichier final : databases/outgroups_all_nr.fa
# ================================================================
# DÉPENDANCES :
#   - vsearch  : déréplication et détection de chimères
#   - zcat     : décompression des archives .gz
#   - awk      : filtrage et manipulation de texte
#   - python3  : suppression des sous-séquences incluses
# ================================================================
# STRUCTURE DES SORTIES :
#   databases/
#   ├── Lactococcus_nr.fa            ← Cible finale (alignement → étape 2)
#   ├── Lactococcus_raw_extracted.fa ← Archive intermédiaire (avant chimères)
#   ├── Lactococcus_nr_chimeras.fa   ← Séquences rejetées (chimères)
#   ├── outgroups_all_nr.fa          ← Tous les outgroups fusionnés
#   ├── Lactobacillus_nr.fa          ← Outgroup individuel
#   ├── Streptococcus_nr.fa          ← Outgroup individuel
#   └── ... (un fichier par genre outgroup)
# =============================================================

set -euo pipefail

# ─────────────────────────────────────────────
# CHEMINS & GENRES
# ─────────────────────────────────────────────

REFSEQ="$HOME/Tenebrion/data/16S/RefSeq/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz" # Chemin vers BDD RefSeq
SILVA="$HOME/Tenebrion/data/16S/SILVA/silva_nr99_v138.2_toSpecies_trainset.fa.gz" # Chemin vers BDD RefSeq
OUTDIR="$HOME/primer_lactococcus/databases"                                       # Dossier de sortie

TARGET="Lactococcus"                                                                     # Genre cible pour les amorces
OUTGROUPS=(Lactobacillus Streptococcus Enterococcus Pediococcus Weissella Leuconostoc)   # Outgroup de genre proche 

# ─────────────────────────────────────────────
# FONCTIONS UTILITAIRES
# ─────────────────────────────────────────────
log()        { echo "[$(date +%H:%M:%S)] $*"; }           # Affiche l'heure
die()        { echo "ERREUR : $*" >&2; exit 1; }          # Message d'erreur si le script plante
count_seqs() { grep -c "^>" "$1" 2>/dev/null || true; }   # Compte le nombre de séquence dans mes fichiers 

# ─────────────────────────────────────────────
# VÉRIFICATION DES DÉPENDANCES ET FICHIERS
# ─────────────────────────────────────────────
check_deps() {
    log "Vérification des dépendances..."

    # Vérifie que chaque outil nécessaire est accessible dans le PATH
    for cmd in vsearch zcat awk python3; do
        command -v "$cmd" >/dev/null 2>&1 \
            || die "$cmd introuvable — vérifie que l'environnement conda est activé."
    done

    # Vérifie que les fichiers sources existent avant de commencer
    [[ -f "$REFSEQ" ]] || die "RefSeq introuvable : $REFSEQ"
    [[ -f "$SILVA"  ]] || die "SILVA introuvable  : $SILVA"
    log "  -> OK"
}
# ─────────────────────────────────────────────
# FONCTION PRINCIPALE : EXTRACTION ET DÉRÉPLICATION
# ─────────────────────────────────────────────
# Arguments :
#   $1 = nom du genre (ex: "Lactococcus", "Streptococcus")
#   $2 = longueur minimale des séquences à conserver (en pb)
#
# Étapes internes :
# - Extraction des séquences du genre depuis RefSeq + SILVA
# - Fusion des deux sources dans un fichier temporaire
# - Filtre ATCG strict : supprime toute séquence contenant des bases ambiguës (N, K, Y, W, R…)
# - Déréplication 100% + filtre longueur minimale (vsearch)
# ─────────────────────────────────────────────
extract_and_derep() {
    local genus="$1"
    local minlen="$2"

    # Répertoire temporaire isolé pour ne pas mélanger les genres
    local tmpdir="$OUTDIR/tmp_${genus}"

    # Fichier de sortie final pour ce genre
    local outfile="$OUTDIR/${genus}_nr.fa"
    mkdir -p "$tmpdir"

    log "  Extraction : $genus (minlen=${minlen} bp)"

    ## Extraction RefSeq
    # zcat décompresse le fichier .gz en flux
    # awk filtre les séquences dont l'en-tête contient le nom du genre :
    #   - Dès qu'une ligne ">" est vue, on remet found à 0
    #   - Si la ligne contient le genre, on active found
    #   - On imprime les lignes où found = 1 (en-tête + séquence)
    # sed ajoute un préfixe "Genre|RefSeq|" à chaque en-tête pour tracer l'origine
    zcat "$REFSEQ" | awk -v g="$genus" '
        /^>/ { found = 0 }
        $0 ~ g { found = 1 }
        found  { print }
    ' | sed "s/^>/>${genus}|RefSeq|/" > "$tmpdir/${genus}_refseq.fa"

    ## Extraction SILVA (pareil que pour RefSeq)
    zcat "$SILVA" | awk -v g="$genus" '
        /^>/ { found = 0 }
        $0 ~ g { found = 1 }
        found  { print }
    ' | sed "s/^>/>${genus}|SILVA|/" > "$tmpdir/${genus}_silva.fa"

    local n_ref=$(count_seqs "$tmpdir/${genus}_refseq.fa")
    local n_sil=$(count_seqs "$tmpdir/${genus}_silva.fa")

    # Fusion des deux BDDS
    cat "$tmpdir/${genus}_refseq.fa" "$tmpdir/${genus}_silva.fa" > "$tmpdir/${genus}_merged.fa"
    local n_merged=$(( n_ref + n_sil ))

    # Filtre ATCG strict
    # Utilise awk pour ne conserver que les séquences constituées
    # UNIQUEMENT des bases A, T, C, G (majuscules ou minuscules).
    # Toute séquence contenant N, K, Y, W, R ou autre code IUPAC
    # est supprimée, car ces ambiguïtés fausseraient le calcul de Tm
    # et la détection de mismatches lors du criblage de spécificité.
    # 
    # Fonctionnement du bloc awk :
    #   RS=">" : chaque enregistrement est délimité par ">"
    #   FS="\n": les champs sont séparés par des sauts de ligne
    #   NR>1   : on saute le premier enregistrement (vide avant le 1er ">")
    #   seq    : reconstruit la séquence en concaténant les lignes 2..NF
    #   Le test regex /^[ATCGatcg]+$/ garantit 100% de bases valides
    awk 'BEGIN {RS=">"; FS="\n"} NR>1 {seq=""; for(i=2;i<=NF;i++) seq=seq $i; if(seq ~ /^[ATCGatcg]+$/) printf ">%s", $0}' \
        "$tmpdir/${genus}_merged.fa" > "$tmpdir/${genus}_atcg_only.fa"
    
    local n_atcg=$(count_seqs "$tmpdir/${genus}_atcg_only.fa")

    ## Déréplication 100% identité + filtre longueur
    # --derep_fulllength : regroupe les séquences identiques à 100% en une seule copie
    # --minseqlength     : élimine les séquences trop courtes (1200 pour le groupe cible et 900 pour les outgroups)
    # --fasta_width 0    : désactive le retour à la ligne automatiquement (une séquence = une seule ligne dans le FASTA)
    vsearch \
        --derep_fulllength "$tmpdir/${genus}_atcg_only.fa" \
        --output           "$outfile" \
        --fasta_width      0 \
        --minseqlength     "$minlen" \
        --quiet

    local n_nr=$(count_seqs "$outfile")

    # Affiche un résumé tabulaire du nettoyage pour ce genre
    printf "    %-15s | Brutes: %4d | ATCG-only: %4d | Uniques (>=%dbp): %4d\n" \
        "$genus" "$n_merged" "$n_atcg" "$minlen" "$n_nr"

    # Nettoyage du répertoire temporaire
    rm -rf "$tmpdir"
}

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
check_deps
mkdir -p "$OUTDIR"
log "Dossier de sortie : $OUTDIR"
echo ""

# =============================================================
# PHASE 1 : EXTRACTION CIBLE ET OUTGROUPS
# =============================================================

## Cible : Lactococcus 
# Seuil de 1200 pb pour garantir des séquences 16S complètes (le gène 16S complet mesure ~1550 pb)
log "=== Constitution de la base cible ==="
extract_and_derep "$TARGET" 1200
echo ""

## Outgroups : 6 genres phylogénétiquement proches
# Seuil de 900 pb accepté pour maximiser la couverture (les séquences partielles sont utiles pour tester la spécificité)
log "=== Constitution des bases outgroups ==="
for genus in "${OUTGROUPS[@]}"; do
    extract_and_derep "$genus" 900
done
echo ""

## Fusion de tous les outgroups en un seul fichier
log "=== Fusion des outgroups ==="
> "$OUTDIR/outgroups_all_nr.fa"  # Crée ou vide le fichier de destination
for genus in "${OUTGROUPS[@]}"; do
    cat "$OUTDIR/${genus}_nr.fa" >> "$OUTDIR/outgroups_all_nr.fa"
done
log "  -> outgroups_all_nr.fa : $(count_seqs "$OUTDIR/outgroups_all_nr.fa") seq. au total"
echo ""

# =============================================================
# PHASE 2 : PRÉTRAITEMENT AVANCÉ DE LA CIBLE (Lactococcus)
# Cette phase s'applique UNIQUEMENT à la cible
# =============================================================
log "=== Pré-traitement de ($TARGET) ==="

# Nombre de séquences avant le prétraitement avancé
N_BEFORE=$(count_seqs "$OUTDIR/${TARGET}_nr.fa")

# On renomme le fichier extrait pour libérer le nom final "Lactococcus_nr.fa"
# Ce fichier _raw_extracted est gardé comme archive de sécurité
mv "$OUTDIR/${TARGET}_nr.fa" "$OUTDIR/${TARGET}_raw_extracted.fa"

# # ÉTAPE A : Suppression des sous-séquences incluses
# Problème : certaines séquences sont de simples fragments d'une autre.
#   Ex: ATCGATCG est contenu dans GGATCGATCGCC → redondance inutile
#   qui gonfle l'alignement et crée du bruit.
#
# Solution en deux temps :
# - vsearch --derep_fulllength avec --strand both :
#   Déréplique aussi sur le brin complémentaire inverse,
#   éliminant les doublons dont un serait "retourné" par erreur
# - Script Python : détecte si une séquence est entièrement
#   contenue dans une autre (test d'appartenance de chaîne),
#   en traitant les séquences de la plus longue à la plus courte
log "  -> Étape A : Filtrage des sous séquences incluses (vsearch + Python)"
vsearch \
    --derep_fulllength "$OUTDIR/${TARGET}_raw_extracted.fa" \
    --output           "$OUTDIR/${TARGET}_derep_both.fa" \
    --fasta_width      0 \
    --strand           both \
    --quiet

# Les variables d'environnement permettent de les passer dans un environnement python
export FASTA_IN="$OUTDIR/${TARGET}_derep_both.fa"
export FASTA_OUT="$OUTDIR/${TARGET}_nosub.fa"

python3 << 'INNERPY'
import os

fasta_in = os.environ.get("FASTA_IN")
fasta_out = os.environ.get("FASTA_OUT")

# ── Lecture du FASTA en mémoire dans un dictionnaire {id: séquence} ──
seqs = {}
curr_id = None
curr_seq = []

with open(fasta_in, "r") as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if curr_id: 
                seqs[curr_id] = "".join(curr_seq)
            curr_id = line[1:].split()[0]  # Prend le premier mot après ">"
            curr_seq = []
        else:
            curr_seq.append(line.upper())
    if curr_id: 
        seqs[curr_id] = "".join(curr_seq)

# ── Algorithme de suppression des sous-séquences ──
# On trie les séquences par longueur décroissante pour comparer
# en priorité les plus longues (elles ne peuvent pas être des sous-seqs
# de séquences plus courtes → gain de temps)
kept = {}
seq_list = sorted(seqs.items(), key=lambda x: len(x[1]), reverse=True)

for sid, seq in seq_list:
    # Vérifie si cette séquence est entièrement contenue dans
    # une séquence déjà acceptée (test d'appartenance Python 'in')
    is_sub = any(seq in kept_seq for kept_seq in kept.values())
    if not is_sub:
        kept[sid] = seq  # Séquence unique : on la conserve

# ── Écriture du FASTA filtré (60 caractères par ligne) ──
with open(fasta_out, "w") as out:
    for sid, seq in kept.items():
        out.write(f">{sid}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")
INNERPY

N_NOSUB=$(count_seqs "$OUTDIR/${TARGET}_nosub.fa")
log "     $(($N_BEFORE - $N_NOSUB)) sous-séquences éliminées."

## ÉTAPE B : Détection et suppression des chimères
# Une séquence chimérique est une séquence "hybride" artificiellement
# assemblée lors du PCR ou du séquençage à partir de deux séquences
# biologiquement distinctes. Elle n'existe pas dans la nature et
# pourrait fausser l'alignement et le design d'amorces.
#
# --uchime_denovo : algorithme UCHIME en mode "de novo"
#   (pas besoin d'une base de référence, les séquences se servent
#    mutuellement de référence en cherchant les "points de croisement")
# --nonchimeras   : séquences validées → fichier final utilisé en aval
# --chimeras      : séquences rejetées → conservées pour inspection
# --uchimealns    : alignement des chimères détectées (débogage)
log "  -> Étape B : Détection des chimères (vsearch --uchime_denovo)"
touch "$OUTDIR/${TARGET}_nr.fa" "$OUTDIR/${TARGET}_nr_chimeras.fa"

vsearch \
    --uchime_denovo  "$OUTDIR/${TARGET}_nosub.fa" \
    --nonchimeras    "$OUTDIR/${TARGET}_nr.fa" \
    --chimeras       "$OUTDIR/${TARGET}_nr_chimeras.fa" \
    --uchimealns     "$OUTDIR/${TARGET}_nr_chimeras_aln.txt" \
    --fasta_width    0 \
    2> "$OUTDIR/uchime.log" || true
    # "|| true" : on ne bloque pas le script si vsearch retourne un warning
    # (comportement attendu quand aucune chimère n'est détectée)

N_FINAL=$(count_seqs "$OUTDIR/${TARGET}_nr.fa")
N_CHIM=$(count_seqs "$OUTDIR/${TARGET}_nr_chimeras.fa" 2>/dev/null || echo 0)

# Suppression des fichiers intermédiaires devenus inutiles
rm -f "$OUTDIR/${TARGET}_derep_both.fa" "$OUTDIR/${TARGET}_nosub.fa"

# ─────────────────────────────────────────────
# RÉSUMÉ FINAL
# ─────────────────────────────────────────────
echo ""
log "================================================================"
log "RÉSUMÉ GÉNÉRAL DU PIPELINE :"
echo "----------------------------------------------------------------"
log "Bases Outgroups prêtes pour spécificité :"
printf "  %-25s %4d séquences\n" "outgroups_all_nr.fa" "$(count_seqs "$OUTDIR/outgroups_all_nr.fa")"
echo ""
log "Cible ($TARGET) prête pour alignement :"
printf "  %-40s %4d séquences\n" "1. Entrée post-extraction"         "${N_BEFORE}"
printf "  %-40s %4d séquences\n" "2. Après filtre sous-séquences"    "${N_NOSUB}"
printf "  %-40s %4d séquences\n" "3. Chimères supprimées"            "${N_CHIM}"
printf "  %-40s %4d séquences\n" "-> FICHIER FINAL (${TARGET}_nr.fa)" "${N_FINAL}"
echo "================================================================"