import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy import constants

from Potentials import pot, discretiser_period, pot_solve
from Solvers import find_eigenvalues, solve, solve_discretized



def main(plots = False): 
    """
    Fonction principale qui calcule et affiche le potentiel périodique d'un cristal 1D.
    
    Cette fonction:
    1. Définit les constantes physiques du système
    2. Discrétise le potentiel périodique
    3. Trace et compare le potentiel continu et sa version discrétisée
    
    Constantes physiques:
    - L : Période du réseau [m]
    - N : Nombre de puits
    - a : Largeur caractéristique des puits [m]
    - e : Charge élémentaire [C]
    - eps0 : Permittivité du vide [F/m]
    """

    #Définition des constante physique :
    L = 5e-10 #m
    a = 1e-10#m
    N= 1
    n = 25
    eps0 = constants.epsilon_0
    e = constants.elementary_charge
    Q = e
    params = (Q, e, a, eps0) 
    x_min = 0.0
    x_max = (N+1)*L
    x = np.linspace(x_min,x_max,1000)

    #Discretisation du potentiel 
    x_bins, V_bins, x_centers = discretiser_period(x_min,x_max,n,pot,L,N,params)
    print(V_bins, x_bins)
    E_min = np.min(V_bins)
    E_max = 0
    E_ok = find_eigenvalues(V_bins, x_bins , E_min, E_max,N,L, plot= True)
    #psi = solve_discretized(V_bins, x_bins, E_ok[0], x, N, L, plot = True)
    #norm = np.trapz((psi)**2, x)
    #if abs(norm - 1.0) < 0.05:  # Tolérance de 5%
    #    print("Test de normalisation RÉUSSI ✓")
    #else:
    #       print("Test de normalisation ÉCHOUÉ ✗")
    
    if plots == True:
        plt.figure(figsize=(10,8))
        plt.plot(x, [pot(xi,params,N,L) for xi in x], color = 'b', label='Potential')
        plt.step(x_bins, np.append(V_bins, V_bins[-1]), where='post', color='black', label='Discretized Potential')
        num_energies_to_plot = len(E_ok)
        #for i in range(num_energies_to_plot):
         #   plt.axhline(y=E_ok[i], color='g', linestyle='--')

        plt.xlabel('x (m)')
        plt.ylabel('Potential Energy (J)')
        plt.title('Potential and Energy Levels')
        plt.legend()
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    main(True)
