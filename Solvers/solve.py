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

def plane_wave(x, k):

    return np.exp(1j * k * x)


def OdSchrodinger(x,y,E,V,params,N,L):
    psi, dpsi_dx = y
    d2psi_dx2 = (V(x,params,N,L)-E)*psi*(2*m/(hbar**2))
    return ([dpsi_dx,d2psi_dx2])

def bc(ya, yb):
    return np.array([ya[0], yb[0]])

def solve(V,E,S0,x_mesh,params,N,L,k):
    y_guess = np.zeros((2,len(x_mesh)))
    x0 = L
    sigma = L/2
    psi_guess = S0 * np.exp(-((x_mesh - x0)**2) / (2 * sigma**2)) * np.cos(k * x_mesh)
    dpsi_guess = (-(x_mesh - x0)/sigma**2 * psi_guess 
                  - k * S0 * np.exp(-((x_mesh - x0)**2)/(2 * sigma**2)) * np.sin(k * x_mesh))
    y_guess[0, :] = psi_guess  # ψ(x)
    y_guess[1, :] = dpsi_guess  # ψ'(x)
    sol =solve_bvp(
    fun=lambda x, y: OdSchrodinger(x, y, E, V,params,N,L),
    bc=lambda ya, yb: bc(ya, yb),
    x=x_mesh,
    y=y_guess,
    max_nodes=10000 
)
    return sol


