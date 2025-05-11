from Potentials import discretiser, pot, pot_period, discretiser_period, potential_voisin
from Solvers import solve, build_transfer_matrix, find_eigenvalues, k_from_energy
import numpy as np
import matplotlib.pyplot as plt


def main(): 
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
    # Définition des constantes physiques du système
    L = 0.5e-9        
    N = 1
    x_min = 0.0       
    x_max = (N+1)*L   
    a = 0.1e-9        
    e = 1.60217663e-19    
    Q = e             
    eps0 = 8.8541878188e-12    
    x0 = L            

    # Configuration des paramètres du potentiel
    params = (Q, e, a, x0, eps0)    
    num_intervalles = 101

    # Définition des points pour la visualisation
    x = np.linspace(x_min, x_max, 1000)    

    # Configuration des paramètres énergétiques
    E_min = -10*7 * 1.6e-19     # -10 eV
    E_max = -1e-20              # ≈0 eV
    E = -10*6*1.6e-19          # Energie test
    S0 = 1e-3                   # Amplitude initiale

    # Discrétisation du potentiel périodique
    x_bins, V_bins, x_centers = discretiser_period(x_min, x_max, num_intervalles, pot_period, L, N, params)
    E_bins = np.linspace(E_min, E_max, num_intervalles)
    dx = x_bins[1] - x_bins[0]  # Largeur de chaque intervalle
    E_ok = find_eigenvalues(V_bins, dx, E_min, E_max, num_points=1000, tol=1e-3, plot=False)
    print(f"Valeurs propres trouvées : {E_ok}")
    for i, energy in enumerate(E_ok):
        psi = solve(pot_period, energy, x_bins, params, N, L, plot=True, 
                    save_path=f'Figures/wavefunction_{i}.pdf')
                
    #psi = solve(pot_period, E_ok, x_bins, params, N, L, plot=True,
               # save_path=f'Figures/wavefunction.pdf')
        

    # Position des puits quantiques
    x_wells = x0 + np.arange(N) * L

    # Visualisation du potentiel
    plt.figure(figsize=(16, 10))  
    plt.plot(x,[pot_period(xi, params,N,L) for xi in x], 
            label="Potentiel périodique", color='blue', linewidth=1.5)
    plt.step(x_bins, np.append(V_bins[0], V_bins), 
            where='pre', label="Découpe", color='black', linewidth=1.5)
    plt.scatter(x_centers, V_bins, 
               color='red', s=15, zorder=3, label="Valeur moyenne de V(x) sur une marche")
    
    # Configuration de l'affichage
    plt.xlabel("x [m]", fontsize=14)
    plt.ylabel("V(x) [J]", fontsize=14)
    plt.legend(loc="best", fontsize=12)
    plt.tight_layout()
    
    # Sauvegarde et affichage du résultat
    plt.savefig('Figures/potentiel.pdf', dpi=300, bbox_inches='tight')
    plt.show()  

if __name__ == "__main__":
    main()

