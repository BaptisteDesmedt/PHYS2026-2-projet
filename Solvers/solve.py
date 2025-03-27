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

def OdSchrodinger(x,y,E,V,params):
    psi, dpsi_dx = y
    d2psi_dx2 = (V(x,params)-E)*psi*(2*m/(hbar**2))
    return ([dpsi_dx,d2psi_dx2])

def bc(ya, yb):
    return np.array([ya[0], yb[0]])

def solve(V,E,S0,x_mesh,params):
    y_guess = np.zeros((2,len(x_mesh)))
    sol =solve_bvp(
    fun=lambda x, y: OdSchrodinger(x, y, E, V,params),
    bc=lambda ya, yb: bc(ya, yb),
    x=x_mesh,
    y=y_guess,
)
    return sol


