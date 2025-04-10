from Potentials import discretiser, pot, pot_period, discretiser_period, potential_voisin
from Solvers import solve, build_transfer_matrix, find_eigenvalues
import numpy as np
import matplotlib.pyplot as plt


def main(): 
    #constante :
    L = 0.5e-9
    N = 2
    x_min = 0.0
    x_max = L*N+ L
    a = 0.1e-9
    e = 1.60217663e-19
    Q = e
    eps0 = 8.8541878188e-12
    x0 = L
    params = (Q, e, a,x0, eps0)
    num_intervalles = 101
    x = np.linspace(x_min, x_max, 1000)
    E_min = -10*7 * 1.6e-19     # -10 eV
    E_max = -1e-20
    E = -10*6*1.6e-19
    S0 = 1e-3
    #x_bins, V_bins, x_centers = discretiser_period(x_min, x_max, num_intervalles,pot_period,L, N, params)
    k = np.pi / L 

    #sol = solve(pot_period, E, S0,x,params,N,L,k)
    x_wells = x0 + np.arange(N) * L

    
    plt.figure(figsize=(14, 6))
    #plt.plot(x,[pot_period(xi, params,N,L) for xi in x], label="Potentiel", color="green", linewidth=1)
    plt.plot(x,[potential_voisin(xi, params,N,L) for xi in x], label="Potentiel 1er voisin", color="blue", linewidth=1)
    #plt.step(x_bins, np.append(V_bins[0], V_bins), where='pre', label="Découpe", color="blue", linewidth=1)
    plt.xlabel("x [m]")
    plt.ylabel("V(x) [J]")
    #plt.plot(x, sol.sol(x)[0], label='Solution numérique')
    plt.xlabel("x [m]")
    plt.ylabel("psi(x)**2 [m]")
    plt.legend(loc = "best")
    plt.show()  
if __name__ == "__main__":
    main()

