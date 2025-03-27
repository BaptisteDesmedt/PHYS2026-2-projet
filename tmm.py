import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar
import pot
import discret

#CONSTANTS DEFINITION

hbar = 1.0545718e-34  # Reduced Planck's constant in J·s
m = 9.10938356e-31    # Electron mass in kg (We assume we are working with electrons)


# Get the wave vector k from the energy E and the potential V
def k_from_energy(E, V):
    if E >= V:
        return np.sqrt(2 * m * (E - V)) / hbar
    else:
        # Return the imaginary wave vector
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

def boudary_cond(V,E,dx):
    M = build_transfer_matrix(V, E, dx)
    return M[0,0]       


def find_eigenvalues(V, dx, E_min, E_max, num_points=1000, tol=1e-3, plot=False):
    # Coarse grid search
    E_values = np.linspace(E_min, E_max, num_points)
    transmission = np.zeros_like(E_values, dtype=complex)
    
    # Store potential eigenvalues
    potential_eigenvalues = []
    
    for i, E in enumerate(E_values):
        transmission[i] = boudary_cond(V, E, dx)
        if np.abs(transmission[i]) < tol:
            potential_eigenvalues.append(E)
    
    # If we found eigenvalues, return them
    if potential_eigenvalues:
        return potential_eigenvalues
    
    # If no eigenvalues found with the tolerance, find the closest ones
    abs_transmission = np.abs(transmission)
    min_indices = np.argsort(abs_transmission)[:5]  # Get indices of 5 smallest values
    
    refined_eigenvalues = []
    for idx in min_indices:
        E_guess = E_values[idx]
        try:
            result = root_scalar(
                lambda E: np.abs(boudary_cond(V, E, dx)),
                bracket=[max(E_min, E_guess-0.1*(E_max-E_min)), min(E_max, E_guess+0.1*(E_max-E_min))],
                method='brentq'
            )
            if result.converged and np.abs(result.root) < 1e-2:
                refined_eigenvalues.append(result.root)
        except ValueError:
            continue
    
    if refined_eigenvalues:
        return refined_eigenvalues
    
    # If still no eigenvalues found, return the energies with minimum transmission
    return E_values[min_indices]
        
    

    
