import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar
from Potentials.pot import pot
from Potentials.discret import discretiser

#CONSTANTS DEFINITION

hbar = 1.0545718e-34  # Reduced Planck's constant in J·s
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
    return M[0,0]       


def find_eigenvalues(V, dx, E_min, E_max, num_points=1000, tol=1e-3, plot=False):
    """
    Finds the eigenvalues of the system within the given energy range.

    Parameters:
    V : array-like
        The potential profile.
    dx : float
        The spatial step size.
    E_min : float
        The minimum energy to search for eigenvalues.
    E_max : float
        The maximum energy to search for eigenvalues.
    num_points : int, optional
        The number of points in the coarse grid search (default is 1000).
    tol : float, optional
        The tolerance for identifying eigenvalues (default is 1e-3).
    plot : bool, optional
        Whether to plot the transmission (default is False).

    Returns:
    list
        A list of eigenvalues. If no eigenvalues are found within the tolerance,
        the function returns the energies with the minimum transmission as a fallback.
    """
    # Create a coarse grid of energies
    E_values = np.linspace(E_min, E_max, num_points)
    transmissions = np.zeros(num_points)

    # Calculate the transmission for each energy value
    for i, E in enumerate(E_values):
        T = boundary_cond(V, E, dx)
        transmissions[i] = np.abs(T)
    
    # Find minima in transmission function (candidates for eigenvalues)
    potential_eigenvalues = []
    mins, _ = find_peaks(-transmissions)
    
    # For each minimum, refine the search to get precise eigenvalue
    if len(mins) > 0:
        for min_idx in mins:
            # Only consider significant minima
            if transmissions[min_idx] < tol:
                # Define a narrower range around the minimum
                E_search_min = E_values[max(0, min_idx-5)]
                E_search_max = E_values[min(num_points-1, min_idx+5)]
                
                # Use root-finding to get precise eigenvalue
                try:
                    result = root_scalar(
                        lambda E: np.abs(boundary_cond(V, E, dx)),
                        bracket=[E_search_min, E_search_max],
                        method='brentq'
                    )
                    if result.converged:
                        potential_eigenvalues.append(result.root)
                except:
                    # Fallback if root finding fails
                    potential_eigenvalues.append(E_values[min_idx])
    
    # If no eigenvalues found, return the energies with minimum transmission
    if len(potential_eigenvalues) == 0 and len(mins) > 0:
        potential_eigenvalues = [E_values[idx] for idx in mins]
    
    # Plot the results if requested
    if plot:
        plt.figure(figsize=(10, 6))
        plt.semilogy(E_values, transmissions)
        plt.scatter([E_values[idx] for idx in mins], [transmissions[idx] for idx in mins], c='r', marker='o')
        plt.scatter(potential_eigenvalues, [tol/2]*len(potential_eigenvalues), c='g', marker='x')
        plt.grid(True)
        plt.xlabel('Energy (J)')
        plt.ylabel('|M[0,0]|')
        plt.title('Finding eigenvalues from transmission minima')
        plt.legend(['Transmission', 'Minima', 'Eigenvalues'])
        plt.tight_layout()
        plt.savefig('Figures/eigenvalues.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    return potential_eigenvalues




