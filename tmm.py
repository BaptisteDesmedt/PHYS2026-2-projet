import numpy as np
import matplotlib.pyplot as plt


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
    
    k = np.zeros_like(V)
    T = np.array([[1,0],
                  [0,1]], dtype=complex)
    k_left = k_from_energy(E, V[0])
    for i in range(len(V)):
        k_right = k_from_energy(E, V[i +1])
        
        M = step_matrix(k_left, k_right)
        T = M @ T
        
        M = propagation_matrix(k_right, dx)
        T = M @ T
        
        k_left = k_right
    k_right = k_from_energy(E, V[-1])
    M = step_matrix(k_left, k_right)
    T = M @ T
    return T 


        
        
        
    

    
