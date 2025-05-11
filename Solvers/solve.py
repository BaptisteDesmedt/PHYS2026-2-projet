#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 19:03:24 2025

@author: baptiste
"""

hbar = 1.0545718e-34  # Reduced Planck's constant in J·s
m = 9.10938356e-31    # Electron mass in kg (We assume we are working with electrons)
import numpy as np 
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
from Potentials.pot import pot_period


def OdSchrodinger(x, y, E, V, params, N, L):
    psi, dpsi_dx = y
    d2psi_dx2 = (V(x, params, N, L) - E) * psi * (2 * m / (hbar**2))
    return ([dpsi_dx, d2psi_dx2])

def bc(ya, yb):
    return np.array([ya[0], yb[0]])

def generate_wavefunction_guess(x_mesh, params, N, L):
    """
    Generate an initial guess for the wavefunction that has a peak at L/2 and is normalized.
    
    Parameters:
    -----------
    x_mesh : array_like
        The x coordinates where the solution is evaluated
    params : tuple
        Parameters for the potential function
    N : int
        Number of periods
    L : float
        Period length
    
    Returns:
    --------
    tuple
        (wavefunction, derivative of wavefunction)
    """
    # Create a Gaussian peak centered at L/2
    center = L/2
    width = L/10  # Width of Gaussian
    
    # Generate Gaussian wavefunction
    psi = np.exp(-((x_mesh - center)**2) / (2 * width**2))
    
    # Normalize the wavefunction
    norm = np.sqrt(np.trapz(psi**2, x_mesh))
    psi = psi / norm
    
    # Calculate derivative using gradient
    dpsi = np.gradient(psi, x_mesh)
    
    return psi, dpsi

def solve(V, E, x_mesh, params, N, L, plot=False, save_path=None):
    """
    Solve the Schrodinger equation.
    
    Parameters:
    -----------
    V : callable
        The potential function
    E : float
        The energy eigenvalue
    x_mesh : array_like
        The x coordinates where the solution is evaluated
    params : tuple
        Parameters for the potential function
    N : int
        Number of periods
    L : float
        Period length
    plot : bool, optional
        Whether to plot the solution and initial guess (default is False)
    save_path : str, optional
        Path to save the figure if plot is True (default is None)
        
    Returns:
    --------
    sol : object
        The solution object from solve_bvp
    """
    # Generate initial guess
    psi_guess, dpsi_guess = generate_wavefunction_guess(x_mesh, params, N, L)
    
    y_guess = np.zeros((2, len(x_mesh)))
    y_guess[::] = psi_guess 
    print("Initial guess for the wavefunction:", y_guess)
    print("Initial guess for the derivative of the wavefunction:", dpsi_guess)
    
    print("Solving the Schrödinger equation...")
    sol = solve_bvp(
        fun=lambda x, y: OdSchrodinger(x, y, E, V, params, N, L),
        bc=lambda ya, yb: bc(ya, yb),
        x=x_mesh,
        y=y_guess,
    )
    if sol.success:
        print("The solution was found successfully!")
        print("solution:", sol)
    else:
        print("The solver failed.")
        exit(1)
    
    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(x_mesh, (sol.sol(x_mesh)[0])**2, label='Wavefunction')
        plt.title(f'Wavefunction for E = {E}')
        plt.xlabel('x')
        plt.ylabel('ψ**2(x)')
        plt.legend()
        plt.grid()
        if save_path is not None:
            plt.savefig(save_path)
        plt.show()
    return sol


