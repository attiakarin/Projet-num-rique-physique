import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Déclaration des variables
E, V0, m, hbar, a = sp.symbols('E V0 m hbar a', real=True, positive=True)
A1 = sp.symbols('A1')
i = sp.I

# Définition de k et q
k = sp.sqrt(2 * m * E) / hbar
q = sp.sqrt(2 * m * (E - V0)) / hbar  # Cas E > V0

# Expression de A3
A3 = (sp.exp(-i * k * a) / 2) * (
    A1 * ((1 + k/q) * sp.exp(i * q * a) + (1 - k/q) * sp.exp(-i * q * a)) +
    A1 * ((k**2 - q**2) * (sp.exp(i * q * a) - sp.exp(-i * q * a)) /
          ((k - q)**2 * sp.exp(i * q * a) - (k + q)**2 * sp.exp(-i * q * a))) *
    ((1 - k/q) * sp.exp(i * q * a) + (1 + k/q) * sp.exp(-i * q * a))
)

# Coefficient de transmission T
T = sp.simplify(sp.Abs(A3)**2 / sp.Abs(A1)**2)

# Substitution des constantes
m_val = 1
hbar_val = 1
a_val = 1
V0_val = 10

# Création de la fonction numérique T(E)
T_func = sp.lambdify(E, T.subs({m: m_val, hbar: hbar_val, a: a_val, V0: V0_val}), 'numpy')

# Valeurs d'énergie (E > V0)
E_vals = np.linspace(V0_val + 0.1, V0_val + 20, 500)
T_vals = T_func(E_vals)

# Tracé du graphe
plt.figure(figsize=(10, 6))
plt.plot(E_vals, T_vals, label=r'$T(E) = \left|\frac{A_3}{A_1}\right|^2$')
plt.xlabel("Énergie E")
plt.ylabel("T(E)")
plt.title("Coefficient de transmission pour E > V₀")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()