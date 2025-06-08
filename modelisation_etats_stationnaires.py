import numpy as np
import matplotlib.pyplot as plt

# Constantes physiques réalistes
hbar = 1.0545718e-34  # J.s
m = 9.10938356e-31    # Masse de l'électron en kg
eV = 1.60218e-19     # 1 eV en joules

# Paramètres du puits
V0 = 10 * eV         # Profondeur du puits en joules
a = 1e-10            # Demi-largeur du puits en mètres (1 Ångström)

# Énergie choisie
E = 50.0 * eV        # Énergie de la particule en joules

# Constantes d'onde
k1 = np.sqrt(2 * m * E) / hbar  # Région I (avant le puits)
q2 = np.sqrt(2 * m * (E + V0)) / hbar  # Région II (dans le puits)
k3 = np.sqrt(2 * m * E) / hbar         # Région III (après le puits)

# On pose A1 = 1
A1 = 1

# Calcul de A2 à partir des conditions de continuité de la fonction d'onde
num_A2 = 2j * k1 * np.exp(-1j * k1 * a)   # Numérateur de A2
den_A2 = (1j * q2 * (np.exp(-1j * q2 * a) - np.exp(3j * q2 * a) * (2j * q2 / (1j * k1 + 1j * q2) - 1)) +
          1j * k1 * (np.exp(-1j * q2 * a) + np.exp(3j * q2 * a) * (2j * q2 / (1j * k1 + 1j * q2) - 1)))
A2 = num_A2 / den_A2  # Calcul de A2 en fonction de A1

# Calcul de A3 à partir de la continuité de la fonction d'onde à x = a
A3 = (2j * q2 / (1j * k1 + 1j * q2)) * A2 * np.exp(1j * q2 * a - 1j * k1 * a)

# Modifié : Calcul de B2 (comme dans l'exemple que vous avez donné)
B2 = A2 * np.exp(2j * q2 * a) * (2j * q2 / (1j * k1 + 1j * q2) - 1)

# Calcul de B1 à partir des conditions de continuité à x = -a
B1 = A1 - A2 * np.exp(-2j * k1 * a) - B2 * np.exp(-2j * k1 * a)

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
x = np.linspace(-8 * a, 8 * a, 1000)

# Puits de potentiel
V = np.zeros_like(x)
V[np.abs(x) <= a] = -V0  # Puits de potentiel dans la région -a à a

plt.figure(figsize=(10, 6))

# Affichage du puits de potentiel
plt.plot(x, V/V0, label=r'Puits de potentiel $V(x)$', color='black', linewidth=2)

# Affichage des fonctions d'onde

# Région I (avant le puits)
plt.plot(x[x < -a], phi_1(x[x < -a], A1=1, B1=B1), label=r'$\phi_1(x)$', color='blue')

# Région II (dans le puits)
plt.plot(x[(x >= -a) & (x <= a)], phi_2(x[(x >= -a) & (x <= a)], A2=A2, B2=B2), label=r'$\phi_2(x)$', color='green')

# Région III (après le puits)
plt.plot(x[x > a], phi_3(x[x > a], A3=A3), label=r'$\phi_3(x)$', color='red')

# Paramètres du graphique
plt.axvline(x=-a, color='black', linestyle='--', label='x = -a')
plt.axvline(x=a, color='black', linestyle='--', label='x = a')
plt.title("Fonctions d'onde dans Puits de potentiel carré fini")
plt.xlabel('Position x (m)')
plt.ylabel('Partie réelle de $\phi(x)$ et $V(x)$')
plt.legend(loc='best')
plt.grid(True)

# Affichage du graphique
plt.show()
