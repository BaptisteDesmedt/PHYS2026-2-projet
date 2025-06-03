import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from Potentials.pot import pot
from Potentials.discret import discretiser

#CONSTANTS DEFINITION
hbar= constants.hbar
m = constants.electron_mass# Electron mass in kg 


# Get the wave vector k from the energy E and the potential V
def k_from_energy(E, V):
    if E >= V:
        return np.sqrt(2 * m * (E - V)) / hbar
    else:
        
        return np.sqrt(2 * m * np.abs((E - V))) / hbar
    
def propagation_matrix(k, dx):
    M0 = np.array([[np.exp(1j * k * dx), 0], 
                   [0, np.exp(-1j * k * dx)]], dtype=complex)
    return M0

def step_matrix(k_left, k_right,x_interface):
    delta = k_left/k_right
    Ms = 0.5*np.array([[-1j*np.exp(-1j*k_right*x_interface)*np.exp(1j*k_left*x_interface)*(1 + delta),
                         1j*np.exp(-1j*k_right*x_interface)*np.exp(-1j*k_left*x_interface)*(delta -1)], 
                   [np.exp(1j*k_right*x_interface)*np.exp(1j*k_left*x_interface)*(1 - delta),
                     -np.exp(1j*k_right*x_interface)*np.exp(-1j*k_left*x_interface)*(1 + delta)]], dtype=complex)
    return Ms

def build_transfer_matrix(V, E, dx_values):
    """
    Build a transfer matrix with variable step sizes.
    
    Parameters:
    -----------
    V : ndarray
        Discretized potential values
    E : float
        Energy
    dx_values : ndarray
        Array of distance intervals between adjacent points
    """
    T = np.eye(2, dtype=complex)
    k_left = k_from_energy(E, V[0])
    
    for i in range(len(V) - 1):
        k_right = k_from_energy(E, V[i+1])
        #print("This is k_left :", k_left," and this is k_right :", {k_right})
        # Interface matrix
        M_interface = step_matrix(k_left, k_right,dx_values[i+1])
        np.dot(M_interface, T, out=T)
        
        # Propagation matrix with variable dx
        #M = propagation_matrix(k_right, dx_values[i])
        #np.dot(M, T, out=T)
        
        
        k_left = k_right
        
    # Final interface
    k_right = k_from_energy(E, V[-1])
    M = step_matrix(k_left, k_right,dx_values[-1])
    np.dot(M, T, out=T)
    
    return T

def boundary_cond(V,E,dx,N,L):
    M = build_transfer_matrix(V, E, dx)
    k = k_from_energy(E, V[-1])
    cond = (M[0,0] - M[0,1])*(np.exp(1j*k*(N+1)*L)) + (M[1,0] - M[1,1])*(np.exp(-1j*k*(N+1)*L))
    return np.abs(cond.real)

def find_eigenvalues(V, dx_values, E_min, E_max,N,L, num_points=10000, tol=1e-19, plot=False):
    E_grid = np.linspace(E_min, E_max, num_points)
    bc_values = np.array([boundary_cond(V, E, dx_values,N,L) for E in E_grid])
    eigenvalues = []

    
    for i in range(len(bc_values) - 1):
        if bc_values[i] * bc_values[i+1] < 0:
            # Linear interpolation for better estimate
            E_zero = E_grid[i] - bc_values[i] * (E_grid[i+1] - E_grid[i]) / (bc_values[i+1] - bc_values[i])
            eigenvalues.append(E_zero)

    if plot == True:
        plt.figure(figsize=(20, 12))
        plt.plot(E_grid, bc_values, label="Boundary condition", color ='darkblue')
        for E in eigenvalues:
            plt.axvline(E, color='r', linestyle='--', alpha=0.5)
        plt.xlabel("Energy (J)")
        plt.ylabel("Boundary condition")
        plt.legend()
        plt.grid
        plt.show()

    return np.array(eigenvalues)


