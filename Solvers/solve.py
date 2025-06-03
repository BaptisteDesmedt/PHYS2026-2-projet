#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 19:03:24 2025

@author: baptiste
"""


import numpy as np 
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

h=6.626e-34
hbar=h/(2*np.pi)
m = 9.11e-31 

def generate_wavefunction_guess(x_mesh, N, L):
    """
    Generate an initial guess for the wavefunction: a Gaussian centered at L/2.
    """
     # Target peak position
    center = L  
    # Create a localized Gaussian-like function centered at peak_position
    width = L /3.5
    psi = 0.2e5*np.exp(-((x_mesh - center) ** 2) / (2 * width ** 2)) +0.2e5*np.exp(-((x_mesh - center*2) ** 2) / (2 * width ** 2))
    
    # Force boundary conditions (ensure it's exactly 0 at boundaries)
    psi[0] = 0.0
    psi[-1] = 0.0
    
    # Calculate the derivative analytically
    dpsi_dx = np.gradient(psi,x_mesh)
    
    return np.vstack([psi, dpsi_dx])
    return y_guess

def bc(ya, yb):
    return np.array([ya[0], yb[0]])

def schrodinger_ode(x, y, E,V,params,N,L):
    # y[0] = psi, y[1] = psi'
    return np.array((y[1], (((2*m)/(hbar**2))*(V(x,params,N,L) - E)) * y[0]))

def schrodinger_ode_discretized(x, y, E, V_interp):
    # y[0] = psi, y[1] = psi'
    return np.array((y[1], (((2 * m) / (hbar**2)) * (V_interp(x) - E)) * y[0]))

def solve(V, E, x_mesh, params, N, L,plot=True, save_path=None):
    """
    Solve the time-independent 1D Schrödinger equation using solve_bvp.
    Returns the normalized wavefunction ψ(x) on x_mesh.
    """
    y_guess  = generate_wavefunction_guess(x_mesh,N,L)


    sol = solve_bvp(lambda x, y: schrodinger_ode(x, y, E,V,params,N,L), bc, x_mesh, y_guess, verbose = 2 , max_nodes = 10000)
    if not sol.success:
        print("The solver failed.", sol.message)
        exit(1)
       
    if plot:
        plt.figure(figsize=(10, 6))
        #plt.plot(sol.x, psi, label='|ψ(x)|²', color = 'darkblue')
        plt.title(f'Wavefunction for E = {E:.3e} J')
        plt.plot(x_mesh , y_guess[0], label ="intial guess", color ='r')
        plt.xlabel('x (m)')
        plt.ylabel('|ψ(x)|²')
        plt.legend()
        plt.grid()
        if save_path is not None:
            plt.savefig(save_path)
        plt.show()

    return sol

def solve_discretized(V_bins, x_bins, E, x_mesh, N, L, plot=True, save_path=None):
    """
    Solve the time-independent 1D Schrödinger equation for a discretized potential.
    """
    # Define the piecewise constant potential function
    def V_partition(x):
        # Find the bin index for each x value
        indices = np.searchsorted(x_bins, x, side='right') - 1
        indices = np.clip(indices, 0, len(V_bins) - 1)  # Ensure indices are within bounds
        return V_bins[indices]

    # Schrödinger equation for the discretized potential
    def schrodinger_ode_discretized(x, y):
        return np.array((y[1], ((2 * m) / (hbar**2)) * (V_partition(x) - E) * y[0]))

    # Generate initial guess for the wavefunction
    y_guess = generate_wavefunction_guess(x_mesh, N, L)

    # Solve the Schrödinger equation
    sol = solve_bvp(schrodinger_ode_discretized, bc, x_mesh, y_guess, verbose=2, max_nodes=1000000)
    if not sol.success:
        print("The solver failed.", sol.message)
        exit(1)

    norm = np.trapz(abs(sol.sol(sol.x)[0])**2, x_mesh)
    psi = abs(sol.sol(sol.x)[0])/np.sqrt(norm)
    # Plot the results if requested
    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(x_mesh, psi, label='|ψ(x)|²', color='darkblue')
        plt.plot(x_mesh,y_guess[0])
        plt.title(f'Wavefunction for E = {E:.3e} J (Discretized Potential)')
        plt.xlabel('x (m)')
        plt.ylabel('|ψ(x)|²')
        plt.legend()
        plt.grid()
        if save_path is not None:
            plt.savefig(save_path)
        plt.show()

    return psi


