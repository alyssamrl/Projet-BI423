#Projet BI423: Détection d'oriC - Code

#Bibliothèques
import re
import matplotlib.pyplot as plt
import customtkinter as ctk
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


#Etape 1 : Charger la séquence à partir du fichier FASTA
def lire_sequence_fasta(nom_fichier):
    """Lit une séquence FASTA (ignore la première ligne de description)"""
    sequence = ""
    with open(nom_fichier, 'r') as fichier:
        #lit la première ligne et l'ignore
        fichier.readline()
        ligne = fichier.readline()
        while ligne:
            sequence += ligne.strip()  #enlever les espaces et les sauts de ligne
            ligne = fichier.readline()  #lire la ligne suivante
    return sequence.upper()


#Etape 2 : Rechercher les régions riches en A et T
def calcul_regions_AT(sequence):
    """Renvoie les régions riches en A/T en comparant le nombre de G/C au nombre de A/T dans une région"""
    taille_région = 300
    régions_dintérêt = [] #liste pour stocker les régions sélectionnées
    #parcours de la séquence par tranches de 300 bases
    for i in range(0, len(sequence) - taille_région + 1, taille_région):
        region = sequence[i:i+taille_région] #extraction de la sous-séquence
        cpt_AT = 0 #compteur A/T
        cpt_CG = 0 #compteur C/G
        #comptage des types de bases dans la région
        for base in region:
            if base == 'A' or base == 'T':
                cpt_AT += 1
            elif base == 'C' or base == 'G':
                cpt_CG += 1
        #vérifie si la région est riche en A/T
        if cpt_AT >= cpt_CG:
            régions_dintérêt.append((i, i + taille_région, region)) #ajoute les positions de début et de fin de région + la région

    return régions_dintérêt


#Etape 3 : Détecter les motifs DnaA
def motifs_consensus(regions):
    """Renvoie les motifs consensus trouvé dans la séquence s'il y en a"""
    motifs = [
        r"TTATCCACA", r"TTGTCCACA", r"TTATGCACA", r"TTATCCAGA", r"TTATTCACA",
        r"TTATCCAAA", r"TCATTCACA", r"TCATTGACA", r"TTGTGCACA", r"ATATTTAAA"
    ]
    résultats = []
    for debut, fin, region in regions: #on cherche des motifs dans les séquences riches en AT
        for motif in motifs: 
            for match in re.finditer(motif, region): #extraire et localiser les motifs
                position_globale = debut + match.start()  #calcul de la position dans la séquence globale
                motif_trouvé = match.group() #motif trouvé dans la région
                résultats.append((position_globale, motif_trouvé))

    return résultats


# Étape 4 : Calculer le GC-skew
def GC_skew(seq):
    """Calcule le GC-skew d'une séquence"""
    C = seq.count('C')
    G = seq.count('G')
    
    somme = G+C

    if somme == 0: 
        return 0 #Pour éviter de faire une division par 0
     
    return (G-C) / somme #calcul du GC-skew

def GCskewfenetre (seq, fenetre):
    """Calcule l'asymétrie GC dans une fenêtre"""
    sequences = [ ]
    position = [ ]
    for i in range (0, len(seq) - fenetre + 1 , fenetre): #Pour calculer fenetre par fenetre
        seq_fenetre = seq[i:i+fenetre]
        val = GC_skew(seq_fenetre)
        sequences.append(val) 
        position.append(i + fenetre//2) # "fenetre//2" est utilisé pour placer chaque valeur de GC_skew au milieu de sa fenêtre (pour mieux visualiser sur le graphe)
    return position, sequences



#Etape 5 : Détecter l'OriC basé sur le GC-skew
def detect_oriC(positions, sequences, seuil=0.1):
    """Détecte l'oriC en cherchant un changement de signe dans le GC-skew."""
    oriC_positions = []
    for i in range(1, len(sequences)):
        if abs(sequences[i] - sequences[i-1]) > seuil and sequences[i-1] < 0 and sequences[i] > 0:  # Détection d'un passage de négatif à positif et dépassement d'un seuil
            oriC_positions.append(positions[i])  # Enregistre la position du changement
    return oriC_positions


# ---- Programme principal ----

# Étape 1 : Lire la séquence
seq = lire_sequence_fasta('seqprojet.txt')

# Étape 2 : Détecter les régions riches en AT
regions = calcul_regions_AT(seq)

# Étape 3 : Détecter les motifs DnaA dans toute la séquence (comme la fonction les cherche dans toute la séquence)
motifs = motifs_consensus(regions)
print("Motifs DnaA trouvés :", motifs)

# Étape 4 : Calculer le GC-skew
fenetre = 5000  # Taille de la fenêtre pour le calcul du GC-skew
pos_gc, val_gc = GCskewfenetre(seq, fenetre) 

# Étape 5 : Filtrer les OriC détectés
oriC_pos = detect_oriC(pos_gc, val_gc)

# Partie affichage du graphe
plt.figure(figsize=(10, 5))
plt.plot(pos_gc, val_gc, marker='o')
plt.title("GC-skew avec fenêtre de taille " + str(fenetre))
plt.xlabel("Position")
plt.ylabel("GC-skew")
plt.axhline(0, color='gray', linestyle='--')  # tracer une ligne horizontale à y = 0 pour visualiser où le GC-skew change de signe

# Marquer les positions de l'OriC détectées sur le graphique
for i in oriC_pos:
    plt.axvline(x=i, color='red', linestyle='--')

plt.grid(True)
plt.show()

# Affichage des positions des OriC filtrées
if oriC_pos:
    print(f"OriC détecté à la position(s) : {oriC_pos}")
else:
    print("Aucun OriC détecté.")


# ----- Interface -----

def charger_fichier(context):
    """Charge un fichier et lit la séquence"""
    filepath = filedialog.askopenfilename(filetypes=[("FASTA / TXT", "*.fasta *.txt")])
    context["sequence"] = lire_sequence_fasta(filepath) #lit la séquence
    context["fenetre_resultats"].insert("end", "Fichier chargé avec succès\n\n") #affiche que le fichier a bien été chargé

def analyser(context):
    """Analyse la séquence et affiche les résultats dans le graphe"""
    sequence = context.get("sequence") 
    fenetre_resultats = context["fenetre_resultats"] #récupère la séquence depuis le contexte
    fenetre_resultats = context["fenetre_resultats"]  #récupère la zone d’affichage des résultats 
    fenetre_resultats.delete("1.0", "end") # efface le contenu précédent dans la zone de texte

    regions = calcul_regions_AT(sequence) #recherche des régions riches en A/T
    motifs = motifs_consensus(regions) #recherche des motifs consensus
    pos_gc, val_gc = GCskewfenetre(sequence, 5000) #calcul du GCskew avec une fenêtre de 5000
    oriC = detect_oriC(pos_gc, val_gc) #détection de l'origine de réplication

    fenetre_resultats.insert("end", f"Motifs DnaA trouvés : {len(motifs)}\n") #affichage des motifs trouvés dans la zone de résultats
    for i, motif in motifs:
        fenetre_resultats.insert("end", f"  → {motif} à la position {i}\n")
    fenetre_resultats.insert("end", f"\nOriC détecté à : {oriC if oriC else 'Aucun'}\n") #affiche l'ori trouvé et sa position

    if context.get("canvas"):
        context["canvas"].get_tk_widget().destroy() #supprime ancien graphique s'il y en a un

    fig = Figure(figsize=(8, 3.5), dpi=100) #création de la figure matplotlib
    graphe = fig.add_subplot(111) #ajoute le graphe dans la figure
    graphe.plot(pos_gc, val_gc, marker='o') #création du graphe
    graphe.axhline(0, color='gray', linestyle='--')
    for i in oriC:
        graphe.axvline(x=i, color='red', linestyle='--')
    graphe.set_title("GC-skew avec fenêtre de taille " + str(fenetre)) #titre du graphe
    graphe.set_xlabel("Position") #titre de l'axe x
    graphe.set_ylabel("GC-skew") #titre de l'axe y
    graphe.grid(True)

    canvas = FigureCanvasTkAgg(fig, master=context["fenetre"]) #associe le graphe à la fenêtre tkinter
    canvas.draw() #ajout du graphe
    canvas.get_tk_widget().pack(pady=10, fill="both", expand=True) #affiche le graphe
    context["canvas"] = canvas #enregistre la fenêtre

#thèmes de l'interface
ctk.set_appearance_mode("light")
ctk.set_default_color_theme("green")

def interface():
    context = {"sequence": None, "canvas": None,} #stocke la séquence et la figure (vide au début)

    fenetre = ctk.CTk() #création de la fenêtre principale
    context["fenetre"] = fenetre

    fenetre.title("Recherche OriC chez une bactérie") #titre de l'interface
    fenetre.geometry("900x700") #taille de l'interface

    titre = ctk.CTkLabel(fenetre, text="Recherche OriC chez une bactérie", font=ctk.CTkFont(size=20, weight="bold")) #ajout du titre en haut de la fenêtre
    titre.pack(pady=20)

    #bouton pour charger le fichier 
    bouton_charger = ctk.CTkButton(fenetre, text="Charger fichier", command=lambda: charger_fichier(context))
    bouton_charger.pack(pady=10) #affichage du bouton

    #bouton pour analyser la séquence
    bouton_analyser = ctk.CTkButton(fenetre, text="Lancer l’analyse", command=lambda: analyser(context))
    bouton_analyser.pack(pady=10) #affichage du bouton

    fenetre_resultats = ctk.CTkTextbox(fenetre, height=200) #création zone de texte pour afficher les résultas
    fenetre_resultats.pack(padx=20, pady=10, fill="x") #affichage du bouton
    context["fenetre_resultats"] = fenetre_resultats #stocke la fenètre contenant les résultats

    fenetre.mainloop()

interface()