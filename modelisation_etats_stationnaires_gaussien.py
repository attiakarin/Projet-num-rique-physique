import numpy as np
import matplotlib.pyplot as plt

# Constantes physiques réalistes
hbar = 1.0545718e-34  # J.s
m = 9.10938356e-31    # Masse de l'électron en kg
eV = 1.60218e-19     # 1 eV en joules

# Paramètres du puits gaussien
V0 = 10 * eV         # Profondeur du puits en joules
sigma = 1e-10        # Écart-type du puits gaussien en mètres

# Énergie choisie
E = 50.0 * eV        # Énergie de la particule en joules

# Constantes d'onde
k1 = np.sqrt(2 * m * E) / hbar  # Région I (avant le puits)
q2 = np.sqrt(2 * m * (E + V0)) / hbar  # Région II (dans le puits)
k3 = np.sqrt(2 * m * E) / hbar         # Région III (après le puits)

# On pose A1 = 1
A1 = 1

# Calcul de A2 à partir des conditions de continuité de la fonction d'onde
num_A2 = 2j * k1 * np.exp(-1j * k1 * sigma)   # Numérateur de A2
den_A2 = (1j * q2 * (np.exp(-1j * q2 * sigma) - np.exp(3j * q2 * sigma) * (2j * q2 / (1j * k1 + 1j * q2) - 1)) +
          1j * k1 * (np.exp(-1j * q2 * sigma) + np.exp(3j * q2 * sigma) * (2j * q2 / (1j * k1 + 1j * q2) - 1)))
A2 = num_A2 / den_A2  # Calcul de A2 en fonction de A1

# Calcul de A3 à partir de la continuité de la fonction d'onde à x = sigma
A3 = (2j * q2 / (1j * k1 + 1j * q2)) * A2 * np.exp(1j * q2 * sigma - 1j * k1 * sigma)

# Modifié : Calcul de B2 (comme dans l'exemple que vous avez donné)
B2 = A2 * np.exp(2j * q2 * sigma) * (2j * q2 / (1j * k1 + 1j * q2) - 1)

# Calcul de B1 à partir des conditions de continuité à x = -sigma
B1 = A1 - A2 * np.exp(-2j * k1 * sigma) - B2 * np.exp(-2j * k1 * sigma)

# Affichage des résultats des constantes
print(f"A2 = {A2}")
print(f"A3 = {A3}")
print(f"B1 = {B1}")
print(f"B2 = {B2}")

# Définition des fonctions d'onde pour chaque région (partie réelle avec termes exponentiels)
def phi_1(x, A1=1, B1=0):
    return np.real(A1 * np.exp(1j * k1 * x) + B1 * np.exp(-1j * k1 * x))  # Superposition de termes exponentiels (Région I)

def phi_2(x, A2=A2, B2=0):
    return np.real(A2 * np.exp(1j * q2 * x) + B2 * np.exp(-1j * q2 * x))  # Superposition de termes exponentiels (Région II)

def phi_3(x, A3=A3):
    return np.real(A3 * np.exp(1j * k3 * x))  # Superposition avec un seul terme exponentiel (Région III)

# Tracé des fonctions d'onde pour chaque région
x = np.linspace(-5 * sigma, 5 * sigma, 2000)

# Puits de potentiel gaussien
V = -V0 * np.exp(-x**2 / (2 * sigma**2))

plt.figure(figsize=(10, 6))

# Affichage du puits de potentiel gaussien
plt.plot(x, V/V0, label=r'Puits de potentiel gaussien $V(x)$', color='black', linewidth=2)

# Affichage des fonctions d'onde

# Région I (avant le puits)
plt.plot(x[x < -sigma], phi_1(x[x < -sigma], A1=1, B1=B1), label=r'$\phi_1(x)$', color='blue')

# Région II (dans le puits)
plt.plot(x[(x >= -sigma) & (x <= sigma)], phi_2(x[(x >= -sigma) & (x <= sigma)], A2=A2, B2=B2), label=r'$\phi_2(x)$', color='green')

# Région III (après le puits)
plt.plot(x[x > sigma], phi_3(x[x > sigma], A3=A3), label=r'$\phi_3(x)$', color='red')

# Paramètres du graphique
plt.axvline(x=-sigma, color='black', linestyle='--', label='x = -sigma')
plt.axvline(x=sigma, color='black', linestyle='--', label='x = sigma')
plt.title("Fonctions d'onde dans un Puits de potentiel Gaussien")
plt.xlabel('Position x (m)')
plt.ylabel('Partie réelle de $\phi(x)$ et $V(x)$')
plt.legend(loc='best')
plt.grid(True)

# Affichage du graphique
plt.show()
