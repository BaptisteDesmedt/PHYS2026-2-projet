import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar
from Potentials.pot import pot
from Potentials.discret import discretiser

#CONSTANTS DEFINITION

h=6.62606896e-34
hbar=h/(2*np.pi)
m = 9.10938356e-31    # Electron mass in kg (We assume we are working with electrons)


# Get the wave vector k from the energy E and the potential V
def k_from_energy(E, V):
    if E >= V:
        return np.sqrt(2 * m * (E - V)) / hbar
    else:
        
        return 1j * np.sqrt(2 * m * (V - E)) / hbar
    
def propagation_matrix(k, dx):
    M0 = np.array([[np.exp(1j * k * dx), 0], 
                   [0, np.exp(-1j * k * dx)]], dtype=complex)
    return M0

def step_matrix(k_left, k_right):
    delta = k_left/k_right
    Ms = 0.5*np.array([[1 + delta, 1 - delta], 
                   [1 - delta, 1 + delta]], dtype=complex)
    return Ms

def build_transfer_matrix(V, E, dx):
    
    T = np.eye(2, dtype=complex)
    k_left = k_from_energy(E, V[0])
    for i in range(len(V) -1 ):
        k_right = k_from_energy(E, V[i +1])
        
        M_interface = step_matrix(k_left, k_right)
        T = M_interface @ T
        
        M = propagation_matrix(k_right, dx)
        T = M @ T
        
        k_left = k_right
    k_right = k_from_energy(E, V[-1])
    M = step_matrix(k_left, k_right)
    T = M @ T
    return T 

def boundary_cond(V,E,dx):
    M = build_transfer_matrix(V, E, dx)
    return M[1,0]       


def find_eigenvalues(V, dx, E_min, E_max, num_points=1000, tol=1e-10, plot=True, refine=True, max_iter =100):
    E_grid = np.linspace(E_min, E_max, num_points)
    bc_values = np.array([boundary_cond(V, E, dx) for E in E_grid])
    
    # Find candidate intervals with sign changes
    signs = np.sign(bc_values)
    zero_points = E_grid[np.abs(bc_values) < tol]  # Direct hits
    crossing_points = []
    
    # Detect sign crossings between grid points
    for i in range(len(E_grid)-1):
        if signs[i] == 0 or signs[i+1] == 0:
            continue
        if signs[i] != signs[i+1]:
            crossing_points.append((E_grid[i], E_grid[i+1]))
    
    # Refine candidates using Brent's method
    eigenvalues = []
    if refine:
        for interval in crossing_points:
            try:
                result = root_scalar(
                    lambda E: boundary_cond(V, E, dx),
                    bracket=interval,
                    method='brentq',
                    xtol=tol,
                    maxiter=max_iter
                )
                if result.converged:
                    eigenvalues.append(result.root)
            except ValueError:
                continue
                
    # Combine results and remove duplicates
    eigenvalues = np.unique(np.concatenate([zero_points, eigenvalues]))
    
    # Visualization
    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(E_grid, bc_values, label='Boundary Condition')
        plt.axhline(0, color='k', linestyle='--', alpha=0.5)
        plt.scatter(eigenvalues, np.zeros_like(eigenvalues), 
                   color='r', zorder=5, label='Eigenvalues')
        plt.xlabel("Energy")
        plt.ylabel("Boundary Condition Value")
        plt.title("Eigenvalue Search Results")
        plt.legend()
        plt.grid(True)
        plt.show()
    
    return eigenvalues


