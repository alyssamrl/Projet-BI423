import re
import matplotlib.pyplot as plt
import streamlit as st

# Étape 1 : Charger la séquence à partir du fichier FASTA
def lire_sequence_fasta(file):
    content = file.read().decode() #Pour lire correctement le fichier FASTA
    lignes = content.splitlines()
    sequence = ""
    for ligne in lignes:
        if not ligne.startswith(">"):
            sequence += ligne.strip()
    return sequence


# Étape 2 : Rechercher les régions riches en A et T
def calcul_régions_AT(sequence):
    """Renvoie les régions riches en A/T en comparant le nombre de G/C au nombre de A/T dans une région"""
    taille_région = 300
    régions_dintérêt = []
    
    for i in range(0, len(sequence) - taille_région + 1, taille_région):
        region = sequence[i:i+taille_région]
        cpt_AT = 0
        cpt_CG = 0
        for base in region:
            if base == 'A' or base == 'T':
                cpt_AT += 1
            elif base == 'C' or base == 'G':
                cpt_CG += 1
        if cpt_AT >= cpt_CG:
            régions_dintérêt.append((i, i + taille_région, region))  # récupère les positions de début et de fin de région + la région

    return régions_dintérêt

# Étape 3 : Détecter les motifs DnaA
def motifs_consensus(seq):
    motifs = [
        r"TTATCCACA", r"TTGTCCACA", r"TTATGCACA", r"TTATCCAGA", r"TTATTCACA",
        r"TTATCCAAA", r"TCATTCACA", r"TCATTGACA", r"TTGTGCACA", r"ATATTTAAA"
    ]
    résultats = []
    for motif in motifs:
        for match in re.finditer(motif, seq):  # extraire et localiser les motifs
            position_globale = match.start()  # calcul de la position dans la séquence globale
            motif_trouvé = match.group()  # motif trouvé dans la région
            résultats.append((position_globale, motif_trouvé))

    return résultats

# GC-skew
def GC_skew(seq):
    """ Fonction qui calcule le GC-skew d'une séquence """
    C = seq.count('C')
    G = seq.count('G')
    somme = G + C

    if somme == 0:
        return 0  # Pour éviter de faire une division par 0
     
    return (G - C) / somme

def GCskewfenetre(seq, fenetre):
    """Fonction qui calcule l'asymétrie GC dans une fenêtre"""
    sequences = []
    position = []
    for i in range(0, len(seq) - fenetre + 1, fenetre): 
        seq_fenetre = seq[i:i+fenetre]
        val = GC_skew(seq_fenetre)
        sequences.append(val)
        position.append(i + fenetre // 2)  # Pour placer chaque valeur de GC_skew au milieu de sa fenêtre (pour mieux visualiser sur le graphe)
    return position, sequences

# Fonction pour détecter l'OriC basé sur le GC-skew avec un seuil
def detect_oriC_filtree(positions, sequences, seuil=0.1):
    """Détecte l'oriC avec un seuil de changement net dans le GC-skew"""
    oriC_positions = []
    for i in range(1, len(sequences)):
        if abs(sequences[i] - sequences[i-1]) > seuil and sequences[i-1] < 0 and sequences[i] > 0:
            oriC_positions.append(positions[i])
    return oriC_positions

# --- Interface Streamlit ---
st.set_page_config(page_title="Analyseur d'OriC", layout="wide", page_icon="🧬")
st.title("🧬 Détécteur d'OriC ")

st.sidebar.header("Paramètres utilisés")
fenetre = 5000
seuil = 0.1
taille_region = 300

st.sidebar.markdown(f"- Taille de la fenêtre GC-skew : **{fenetre}**")
st.sidebar.markdown(f"- Seuil de détection OriC : **{seuil}**")
st.sidebar.markdown(f"- Taille des régions AT : **{taille_region}**")


fasta_file = st.file_uploader("Chargez votre fichier FASTA")

if fasta_file:
    seq = lire_sequence_fasta(fasta_file)
    st.success(f"Fichier chargé, longueur de la séquence : {len(seq)} bases.")

    # Motifs DnaA
    motifs = motifs_consensus(seq)
    st.subheader("Motifs DnaA trouvés")
    st.write(f"Nombre de motifs détectés : {len(motifs)}")
    if motifs:
        st.dataframe([{"Position": m[0], "Motif": m[1]} for m in motifs])

    # GC-skew
    pos_gc, val_gc = GCskewfenetre(seq, fenetre)
    oriC_pos = detect_oriC_filtree(pos_gc, val_gc, seuil)

    st.subheader("GC-skew et détection de l'OriC")
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(pos_gc, val_gc, marker='o', label="GC-skew")
    ax.axhline(0, color='gray', linestyle='--')
    for i in oriC_pos:
        ax.axvline(x=i, color='red', linestyle='--', label="OriC détecté" if i == oriC_pos[0] else "")
    ax.set_xlabel("Position")
    ax.set_ylabel("GC-skew")
    ax.set_title("GC-skew de la séquence")
    ax.legend()
    st.pyplot(fig)

    if oriC_pos:
        st.success(f"OriC détecté à la position(s) : {oriC_pos}")
    else:
        st.warning("Aucun OriC détecté.")

    st.markdown("---")
    st.info("Développé par Farah YOUNES et Alyssa Morellon")

else:
    st.info("Veuillez charger un fichier FASTA pour commencer l'analyse.")