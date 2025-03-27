import pot 
import discret
import tmm
import numpy as np
import matplotlib.pyplot as plt

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
    x = np.linspace(x_min, x_max, 1000)
    E = 0.0
    E_min = -2e18     # -10 eV
    E_max = -1e-16
    
    
    #discretisation du potentielle 
    x_steps, V_steps, x_bins, V_bins = discret.discretiser(x_min, x_max, num_intervalles, pot.pot, params)
    L = abs(x_bins[1] - x_bins[0])
    #Construction de la matrice de transfer
    T = tmm.build_transfer_matrix(V_bins, E, L)
    
    #Trouver la valeur de l'energie admise
    E_admitted = tmm.find_eigenvalues(V_bins, L, E_min, E_max)


    #Traçage des graph
    plt.plot(x_steps, V_steps, label="Discrete Steps")
    plt.plot(x, pot.pot(x, params), label="Continuous Function", linestyle='-')
    plt.legend()
    plt.show()  
    print(T)
    return x_steps, V_steps, E_admitted

global x_steps, V_steps,E_admitted
x_steps, V_steps, E_admitted = main()

