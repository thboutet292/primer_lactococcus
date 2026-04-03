#!/usr/bin/env python3
# =============================================================
# plot_silva.py
# Script de visualisation : Rapport final de validation des amorces
#
# Ce script génère un rapport graphique PNG résumant les résultats
# de la validation des amorces contre la base de données SILVA.
# Il est appelé après un BLAST ou une recherche de similarité
# entre les amorces choisies et les séquences SILVA, dont les
# résultats sont fournis sous forme d'un fichier CSV.
# USAGE :
#   python3 plot_silva.py \
#       -i resultats_silva.csv \
#       -f_seq ATCGATCG... -f_tm 58.3 -f_gc 52.1 -f_len 20 \
#       -r_seq TAGCTAGC... -r_tm 57.1 -r_gc 50.0 -r_len 21 \
#       -s 350
#
# ARGUMENTS :
#   -i / --input : fichier CSV des résultats SILVA
#                  (séparateur ";" requis)
#                  Colonnes attendues : "path" (taxonomie) + "organismName"
#   -f_seq       : séquence 5'-> 3' de l'amorce Forward
#   -f_tm        : température de fusion de la Forward (°C)
#   -f_gc        : teneur en GC de la Forward (%)
#   -f_len       : longueur de la Forward (nt)
#   -r_seq       : séquence 5'-> 3' de l'amorce Reverse
#   -r_tm        : température de fusion de la Reverse (°C)
#   -r_gc        : teneur en GC de la Reverse (%)
#   -r_len       : longueur de la Reverse (nt)
#   -s / --size  : taille attendue de l'amplicon (pb)
#
# DÉPENDANCES :
#   pip install pandas plotly kaleido
#
# SORTIE :
#   Rapport_Lactococcus_{taille}bp.png
#   (ex: Rapport_Lactococcus_350bp.png)
# =============================================================

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import argparse
import sys

def main():
    # ─────────────────────────────────────────────
    # PARSING DES ARGUMENTS EN LIGNE DE COMMANDE
    # ─────────────────────────────────────────────
    parser = argparse.ArgumentParser(description="Rapport final de validation des amorces Lactococcus")
    parser.add_argument("-i", "--input", required=True, help="Fichier CSV de SILVA")
    # Propriétés de l'amorce Forward
    parser.add_argument("-f_seq", required=True, help="Séquence 5'→3' de la Forward")
    parser.add_argument("-f_tm",  required=True, help="Tm de la Forward (°C)")
    parser.add_argument("-f_gc",  required=True, help="GC% de la Forward")
    parser.add_argument("-f_len", required=True, help="Longueur de la Forward (nt)")
    # Propriétés de l'amorce Reverse
    parser.add_argument("-r_seq", required=True, help="Séquence 5'→3' de la Reverse")
    parser.add_argument("-r_tm",  required=True, help="Tm de la Reverse (°C)")
    parser.add_argument("-r_gc",  required=True, help="GC% de la Reverse")
    parser.add_argument("-r_len", required=True, help="Longueur de la Reverse (nt)")
    # Paramètre global
    parser.add_argument("-s", "--size", required=True, help="Taille estimée de l'amplicon (pb)")
    args = parser.parse_args()

    # ─────────────────────────────────────────────
    # 1. TRAITEMENT DES DONNÉES SILVA
    # ─────────────────────────────────────────────
    # Lecture du CSV de résultats
    # Séparateur ";"
    try:
        df = pd.read_csv(args.input, sep=";")
    except Exception as e:
        print(f"Erreur lecture : {e}"); sys.exit(1)

    total_hits = len(df)                               # Nombre total de séquences trouvées dans SILVA
    target_genus = "Lactococcus"                       # Genre que nos amorces doivent cibler

    # Extraction du genre réel depuis le champ taxonomique "path"
    # Exemple de "path" : "'Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Lactococcus'"
    # L'avant-dernier segment ([-2]) correspond au genre
    df['Genre_Reel'] = df['path'].apply(
        lambda x: x.strip("'").split(';')[-2] if isinstance(x, str) else "Unknown"
    )
    
    def clean_taxonomy(name):
        """
        Nettoie et simplifie le nom taxonomique pour l'affichage.
        
        Stratégie :
        - Si le nom contient des mots-clés environnementaux/ambigus
          ("uncultured", "bacterium", "metagenome"…), on les regroupe
          sous "Lactococcus sp. (uncultured)" pour ne pas polluer le graphique
          avec des centaines de noms différents
        - Sinon, on conserve uniquement les deux premiers mots de l'espèce
        """
        name_lower = str(name).lower()                                         # nom d'organisme brut (str)

        # Mots-clés indiquant une séquence environnementale ou mal caractérisée
        env_keywords = ["uncultured", "bacterium", "metagenome", "unidentified", "enrichment"]
        if any(k in name_lower for k in env_keywords):
            return f"{target_genus} sp. (uncultured)"

        # Conservation de "Genre espèce" uniquement (2 premiers mots)
        parts = str(name).split()
        return " ".join(parts[:2]) if len(parts) >= 2 else name                # nom nettoyé

    df['Taxon_Display'] = df['organismName'].apply(clean_taxonomy)

    # Calcul du score de spécificité
    # Fraction des hits qui appartiennent au genre cible
    n_target   = len(df[df['Genre_Reel'] == target_genus])
    spec_score = (n_target / total_hits) * 100 if total_hits > 0 else 0

    # Top 10 des taxons les plus fréquents
    # Triés par ordre croissant pour l'affichage horizontal (le plus fréquent en haut)
    species_counts = df['Taxon_Display'].value_counts().head(10).sort_values(ascending=True)

    # ─────────────────────────────────────────────
    # 2. CRÉATION DE LA FIGURE PLOTLY
    # ─────────────────────────────────────────────
    # Tableau des propriétés des amorces
    # Histogramme des taxons SILVA
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.35, 0.65],
        vertical_spacing=0.12,
        specs=[[{"type": "table"}], [{"type": "bar"}]]  # Types de graphiques par panneau
    )

    # 1) Tableau récapitulatif des amorces
    fig.add_trace(
        go.Table(
            header=dict(
                values=["<b>CARACTÉRISTIQUE</b>", "<b>AMORCE FORWARD</b>", "<b>AMORCE REVERSE</b>"],
                fill_color='#1f3b4d',  # Fond sombre pour l'en-tête
                align='left',
                font=dict(color='white', size=12)
            ),
            cells=dict(
                values=[
                    # Colonne 1 : noms des propriétés
                    ["Séquence (5'-3')", "Longueur", "Temp. Fusion (Tm)", "Teneur GC%", "<b>Verdict Genre</b>"],
                    # Colonne 2 : valeurs Forward
                    [args.f_seq, f"{args.f_len} bp", f"{args.f_tm} °C", f"{args.f_gc} %", f"Cible : {target_genus}"],
                    # Colonne 3 : valeurs Reverse + score de spécificité global
                    [args.r_seq, f"{args.r_len} bp", f"{args.r_tm} °C", f"{args.r_gc} %", f"Précision : {spec_score:.1f}%"]
                ],
                fill_color='#f9f9f9',
                align='left',
                font=dict(size=11, family="Arial")
            )
        ),
        row=1, col=1
    )

    # ── Panneau 2 : Histogramme horizontal des taxons SILVA ──
    # Chaque barre représente un taxon (nom d'organisme) et sa longueur
    # indique le nombre de hits dans les résultats SILVA
    fig.add_trace(
        go.Bar(
            x=species_counts.values,    # Nombre de hits (axe X)
            y=species_counts.index,     # Noms des taxons (axe Y)
            orientation='h',            # Histogramme horizontal
            marker=dict(
                color='#3498db',             # Bleu pour les barres
                line=dict(color='#2980b9', width=1)  # Contour légèrement plus foncé
            ),
            text=species_counts.values, # Affiche le nombre de hits sur chaque barre
            textposition='outside',     # Texte à l'extérieur (à droite) des barres
        ),
        row=2, col=1
    )

    # ─────────────────────────────────────────────
    # 3. MISE EN PAGE GLOBALE
    # ─────────────────────────────────────────────
    fig.update_layout(
        title=dict(
            text=f"<b>Analyse de Spécificité Génomique : Genre {target_genus}</b>",
            x=0.5, y=0.98,           # Centré en haut
            xanchor='center',
            font=dict(size=19, color='#1f3b4d')
        ),
        template="simple_white",     # Fond blanc épuré
        height=850, width=1000,      # Dimensions en pixels
        showlegend=False,
        margin=dict(t=110, b=50, l=50, r=50)  # Marges internes
    )

    # Légende de l'axe X du panneau histogramme
    fig.update_xaxes(
        title_text=f"Nombre de hits (Total : {total_hits} séquences)",
        row=2, col=1,
        showgrid=True,
        gridcolor='#eee'
    )
    # Noms des taxons en italique (convention taxonomique)
    fig.update_yaxes(
        tickfont=dict(size=11, style='italic'),
        row=2, col=1
    )

    # ─────────────────────────────────────────────
    # 4. SAUVEGARDE DU RAPPORT EN PNG
    # ─────────────────────────────────────────────
    # Le nom du fichier intègre la taille de l'amplicon pour
    # pouvoir comparer plusieurs paires d'amorces facilement
    # scale=3 : multiplie la résolution (300 dpi effectifs)
    output_file = f"Rapport_Lactococcus_{args.size}bp.png"
    fig.write_image(output_file, scale=3)
    
    print(f"\n[SUCCÈS] Rapport généré : {output_file}")

if __name__ == "__main__":
    main()