import pot 
import discret
import tmm
import numpy as np
import matplotlib.pyplot as plt
import solve

def main(): 
    #constante :
    a = 0.01e-9
    e = 1.6e-19
    Q = e
    eps0 = 8.85e-12
    x0 = 0
    params = (Q, e, a, x0, eps0)
    L = 0.05e-9
    x_min = -L
    x_max = L
    num_intervalles = 10
    N = 10
    x = np.linspace(x_min, x_max, 1000)
    E = 0.0
    E_min = -10 * 1.6e-19     # -10 eV
    E_max = -1e-20
    
    
    #discretisation du potentielle 
    x_steps, V_steps, x_bins, V_bins = discret.discretiser(x_min, x_max, num_intervalles, pot.pot, params)
    dx = abs(x_bins[1] - x_bins[0])
    #Construction de la matrice de transfer
    T = tmm.build_transfer_matrix(V_bins, E, dx)
    
    #Trouver la valeur de l'energie admise
    E_admitted = tmm.find_eigenvalues(V_bins, dx, E_min, E_max)
    
    #Résoudre l'eq
    x_mesh = np.linspace(-L, L, 1000)  # Maillage spatial
    S0 = np.array([0,0])
    psi = solve.solve(pot.pot,E_admitted[0],S0,x_mesh, params)


    #Traçage des graph
    
    plt.plot(psi.x, psi.y[0], label="psi")
    plt.legend()
    plt.show()  
    print(T)
    return psi.x, psi.y ,E_admitted
global x_steps, V_steps,E_admitted
x_steps, V_steps, E_admitted = main()

