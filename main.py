import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

from Potentials import pot, discretiser_period, pot_solve
from Solvers import find_eigenvalues, solve

def test_solver():
    """
    Fonction de test pour vérifier que le solveur d'équations différentielles fonctionne correctement.
    Vérifie la résolution pour des cas simples dont on connaît les solutions analytiques:
    1. Puits infini : E_n = (n²π²ħ²)/(2mL²)
    2. Oscillateur harmonique : E_n = ħω(n+1/2)
    """
    print("Exécution des tests du solveur d'équations différentielles...")
    
    # Test 1: Puits de potentiel infini (particule dans une boîte)
    def infinite_well(x, params, N, L):
        width = L
        if 0 <= x <= width:
            return 0
        else:
            return 1e10  # Une valeur très grande pour simuler un mur infini
    
    # Paramètres pour le puits infini
    L_well = 1e-9  # 1 nm
    N_well = 1
    params_well = ()
    x_min_well = 0
    x_max_well = L_well
    n_well = 1000
    
    # Discrétisation du potentiel de puits infini
    x_bins_well, V_bins_well, _ = discretiser_period(x_min_well, x_max_well, n_well, infinite_well, L_well, N_well, params_well)
    dx_well = x_bins_well[1] - x_bins_well[0]
    
    # Calcul des énergies propres numériques
    E_min_well = 0
    E_max_well = 15 * 1.60218e-19  # 15 eV en Joules
    eigenvalues_well = find_eigenvalues(V_bins_well, dx_well, E_min_well, E_max_well, plot=False)
    
    # Calcul des énergies propres analytiques pour un puits infini
    # E_n = (n²π²ħ²)/(2mL²), n = 1, 2, 3, ...
    hbar = 1.0545718e-34  # J·s
    m = 9.10938356e-31    # kg
    analytical_energies_well = [(n**2 * np.pi**2 * hbar**2)/(2 * m * L_well**2) for n in range(1, 4)]
    
    print("Test du puits infini:")
    # Conversion des résultats en eV pour affichage
    print(f"Énergies numériques trouvées (eV): {eigenvalues_well[:3]/1.60218e-19}")
    print(f"Énergies analytiques attendues (eV): {[E/1.60218e-19 for E in analytical_energies_well]}")
    
    # Calcul de l'erreur relative
    if len(eigenvalues_well) >= 3:
        errors_well = [abs(eigenvalues_well[i] - analytical_energies_well[i])/analytical_energies_well[i] 
                       for i in range(min(3, len(eigenvalues_well)))]
        print(f"Erreurs relatives: {errors_well}")
        if all(error < 0.05 for error in errors_well):  # Tolérance de 5%
            print("Test du puits infini RÉUSSI ✓")
        else:
            print("Test du puits infini ÉCHOUÉ ✗")
    else:
        print("Pas assez de valeurs propres trouvées pour le puits infini")
    
    # Test 2: Test de normalisation
    # Vérifier que les fonctions d'onde sont normalisées
    """
    if len(eigenvalues_well) > 0:
        psi_test = solve(infinite_well, eigenvalues_well[0], x_bins_well, params_well, N_well, L_well, plot=False)
        norm = np.trapz(np.abs(psi_test)**2, x_bins_well)
        print(f"\nTest de normalisation:")
        print(f"Norme de la fonction d'onde = {norm}")
        if abs(norm - 1.0) < 0.05:  # Tolérance de 5%
            print("Test de normalisation RÉUSSI ✓")
        else:
            print("Test de normalisation ÉCHOUÉ ✗")
    """
    print("Tests terminés.")
    return

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
    # Exécution des tests du solveur
    test_solver()
    """
    #Définition des constante physique :
    L = 5e-10 #m
    a = 1e-10#m
    N= 1
    n = 101
    e = -1.6022e-19 #C
    Q = e
    K = 8.987551784e9 #N*m**2/C**2 = à 1/4pieps0
    params = (Q, e, a, K) 
    x_min = 0.0
    x_max = (N+1)*L
    x = np.linspace(x_min,x_max, 10000)
    x_bins, V_bins, x_centers = discretiser_period(x_min,x_max,n,pot,L,N,params)
    
    plt.figure(figsize=(10,8))
    plt.plot(x, [pot(xi,params,N,L) for xi in x], color = 'b')
    
    # Fix the step plot - use where='post' to align values with left edges
    plt.step(x_bins[:-1], V_bins, where='post', color = 'black')
    E_min = 2*np.min(V_bisn)
    E_max = 2*np.max(V_bins)
    dx = x_bins[2] - x_bins[1]
    E_ok = find_eigenvalues(V_bins, dx , E_min, E_max)
    print(E_ok)
    psi = solve(pot_solve,E_ok[0],x_bins, params, N, L, plot = True )
    plt.show()
    """


if __name__ == "__main__":
    main()

