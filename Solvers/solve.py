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
from Potentials.pot import pot_solve


def OdSchrodinger(x, y, E, V, params, N, L):
    psi, dpsi_dx = y
    d2psi_dx2 = (V(x, params, N, L) - E) * psi * (2 * m / (hbar**2))
    return ([dpsi_dx, d2psi_dx2])

def bc(ya, yb):
    return np.array([ya[0], yb[0]])

def generate_wavefunction_guess(x_mesh, params, N, L):
    """
    Generate an initial guess for the wavefunction that:
    - Equals zero at both boundaries
    - Has a peak of amplitude 2 at x = 0.5e-9
    
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
    # Target peak position
    peak_position = 0.5e-9
    
    # Create a localized Gaussian-like function centered at peak_position
    width = L/2  # Adjust width as needed for appropriate spread
    
    # Calculate the normalized distance from the peak position
    # (scaled by domain width to handle different domain sizes)
    normalized_distance = (x_mesh - peak_position)**2 / width**2
    
    # Create a function that's 0 at boundaries and has peak at specified position
    psi = 1.41 * np.exp(-normalized_distance)
    
    # Force boundary conditions (ensure it's exactly 0 at boundaries)
    psi[0] = 0.0
    psi[-1] = 0.0
    
    # Calculate the derivative analytically
    dpsi_dx = -2.0 * (x_mesh - peak_position) / width**2 * psi
    
    return psi, dpsi_dx

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
    psi_guess, dpsi_dx_guess = generate_wavefunction_guess(x_mesh, params, N, L)
    
    y_guess = np.zeros((2, len(x_mesh)))
    y_guess[0] = psi_guess
    y_guess[1] = dpsi_dx_guess

    
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
        plt.plot(x_mesh, (sol.y[0])**2, label='Wavefunction')
        #plt.plot(x_mesh, psi_guess**2)
        #plt.ylim(-4, 8)
        plt.title(f'Wavefunction for E = {E}')
        plt.xlabel('x')
        plt.ylabel('ψ**2(x)')
        plt.legend()
        plt.grid()
        if save_path is not None:
            plt.savefig(save_path)
        plt.show()
    return sol


