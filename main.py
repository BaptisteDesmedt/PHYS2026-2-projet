from Potentials import discretiser, pot, pot_period, discretiser_period
from Solvers import solve, build_transfer_matrix, find_eigenvalues
import numpy as np
import matplotlib.pyplot as plt


def main(): 
    #constante :
    L = 0.5e-9
    N = 1
    x_min = 0.0
    x_max = L*N+ L
    a = 0.1e-9
    e = 1.60217663e-19
    Q = e
    eps0 = 8.8541878188e-12
    x0 = L
    params = (Q, e, a,x0, eps0)
    num_intervalles = 11
    x = np.linspace(x_min, x_max, 1000)
    E_min = -10 * 1.6e-19     # -10 eV
    E_max = -1e-20

    
    plt.figure(figsize=(12, 6))
    plt.plot(x,[pot_period(xi, params,N,L) for xi in x], label="Potentiel", color="green", linewidth=1)
    #plt.step(x_bins, np.append(V_bins[0], V_bins), where='pre', label="Découpe", color="blue", linewidth=1)
    #plt.scatter(x_centers , V_bins , color="red", linewidth=0.5, label="mesh")
    plt.xlabel("x [m]")
    plt.ylabel("V(x) [J]")
    plt.legend()
    plt.show()  

if __name__ == "__main__":
    main()

