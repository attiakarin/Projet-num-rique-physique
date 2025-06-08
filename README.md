# Projet-num-rique-physique

Projet de Physique Moderne – Effet Ramsauer–Townsend

PréIng 2 MI5 – CY Tech 2024–2025

Sujet du projet:

Ce travail s’inscrit dans le cadre du cours de Physique Moderne.
Il s'appuie sur le document officiel « Projet_2425.pdf », décrivant l’effet Ramsauer–Townsend et les étapes à suivre.
Le modèle utilisé est un puits de potentiel 1D à profondeur finie, dans lequel on analyse la diffusion d’une onde.

Objectif du projet: 

Ce projet a pour but d’étudier la transmission quantique d’une particule à travers un puits de potentiel carré fini. 
L’objectif est d’expliquer l’effet Ramsauer–Townsend, un phénomène où la probabilité de diffusion d’une particule devient très faible (voire nulle) à certaines énergies.

Description des scripts utilisés:

courbe_transmission_energie.py

Calcule symboliquement le coefficient de transmission T(E) en fonction de l’énergie.

Utilise la bibliothèque SymPy pour manipuler des expressions analytiques.

Produit un graphe T(E) pour visualiser les zones de transmission maximale.

modelisation_etats_stationnaires.py

Simule les fonctions d’onde phi(x) dans un puits de potentiel.

Utilise des valeurs physiques réalistes 
(énergie en eV, masse de l’électron…).

Affiche le puits de potentiel et les trois régions de l’espace avec la forme des ondes.


Utilisation des fichiers : 

Installer les bibliothèques nécessaires.

pip install numpy matplotlib sympy

Lancer l’un des scripts dans un terminal :

python modelisation_etats_stationnaires.py
ou
python courbe_transmission_energie.py

Bibliothèques utilisées:

numpy: pour les calculs numériques
matplotlib : pour afficher les graphiques
sympy: pour les manipulations symboliques
