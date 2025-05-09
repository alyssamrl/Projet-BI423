## Détection de l’OriC dans une séquence bactérienne

Ce projet Python a pour objectif de détecter les origines de réplication (OriC) dans les génomes bactériens en analysant :

* les régions riches en A/T,
* les motifs consensus DnaA,
* l’asymétrie en GC-skew.

Une interface graphique est également fournie pour charger une séquence, lancer l’analyse et visualiser les résultats.


## Fonctionnalités

* Chargement de séquences ADN au format FASTA (.fasta / .txt)
* Détection de régions riches en A/T
* Recherche de motifs consensus DnaA (séquences cibles de la protéine d’initiation de la réplication)
* Calcul et tracé du GC-skew pour repérer l’OriC
* Détection automatique des changements de signe du GC-skew (indicateurs de l’OriC)
* Interface graphique interactive avec affichage des résultats


## Dépendances

Ce projet utilise les bibliothèques suivantes :

```bash
pip install matplotlib customtkinter
```


## Structure du projet

* `YOUNESMORELLON.py` ou `YOUNESMORELLON.ipynb` : script principal contenant toutes les fonctions d’analyse et l’interface utilisateur
* `sequenceprojet.txt` : exemple de fichier de séquence FASTA utilisé dans notre projet


## Lancer le projet

1. Cloner le dépôt ou copier les fichiers du projet.
2. Lancer le script dans le fichier YOUNESMORELLON.py
3. L’interface graphique s’ouvrira :

   * Cliquez sur "Charger fichier" pour importer la séquence 'sequenceprojet.txt'
   * Cliquez sur "Lancer l’analyse" pour démarrer l’analyse.


## ATTENTION

Lors de l'exécution du code dans un fichier .py, l'interface ne s'affiche pas car matplotlib et tkinter ne fonctionnent pas correctement dans un script Python classique sans environnement interactif. C'est pourquoi nous avons ajouté un fichier .ipynb contenant le même code, afin de pouvoir l'exécuter et visualiser l'interface!!!
