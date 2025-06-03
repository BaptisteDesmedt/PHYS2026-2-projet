import numpy as np 
import math

def coulomb_potential(x, Q, x0, a, epsilon0,e):
    numerator = -Q * e
    denominator = 4 * np.pi * epsilon0 * np.sqrt((x - x0)**2 + a**2)
    return numerator / denominator

def pot(x, params, N, L ):
    """
    Calculate the potential at a point x due to a charge distribution.

    Parameters:
    x : float
        The position where the potential is calculated.
    params : tuple
        A tuple containing the parameters (Q, e, a, x0, eps0):
        Q : float
            Charge magnitude.
        e : float
            Elementary charge.
        a : float
            A constant related to the charge distribution.
        x0 : float
            The position of the charge distribution.
        eps0 : float
            Permittivity of free space.

    Returns:
    float
        The potential at position x.
    """
    Q, e, a, EPS0 = params 
    
    
    if x <= L:
        V = coulomb_potential(x,Q,L,a,EPS0,e)
        if N != 1 :
            V += coulomb_potential(x,Q,2*L,a,EPS0,e)
        return V
    if x >= N*L:
        V = coulomb_potential(x,Q,N*L,a,EPS0,e)
        if N != 1 : 
            V += coulomb_potential(x,Q,(N-1)*L, a, EPS0,e)
        return V
    else : 
        indice = math.floor(x/L)
        V = coulomb_potential(x,Q,indice*L, a , EPS0, e)
        V += coulomb_potential(x,Q,(indice+1)*L,a,EPS0,e)
        return V
        
    return 0

def pot_solve(x, params, N, L):
    """
    Vectorized version of pot() function that calculates the potential at array of points x.
    
    Parameters:
    x : float or ndarray
        Point(s) where the potential is calculated.
    params : tuple
        A tuple containing the parameters (Q, e, a, EPS0)
    N : int
        Number of periods
    L : float
        Period length
    
    Returns:
    float or ndarray
        The potential at x.
    """
    Q, e, a, EPS0 = params  
    x_arr = np.atleast_1d(x)
    V = np.zeros_like(x_arr, dtype=float)
    
    # Vectorized coulomb potential 
    def coulomb_vec(x_vec, x0):
        numerator = -Q * e
        denominator = 4 * np.pi * EPS0 * np.sqrt((x_vec - x0)**2 + a**2)
        return numerator / denominator
    
    # Region where x <= L
    idx1 = (x_arr <= L)
    if np.any(idx1):
        V[idx1] = coulomb_vec(x_arr[idx1], L)
        if N != 1:
            V[idx1] += coulomb_vec(x_arr[idx1], 2*L)
    
    # Region where x >= N*L
    idx2 = (x_arr >= N*L)
    if np.any(idx2):
        V[idx2] = coulomb_vec(x_arr[idx2], N*L)
        if N != 1:
            V[idx2] += coulomb_vec(x_arr[idx2], (N-1)*L)
    
    # Region in between
    idx3 = ~(idx1 | idx2)
    if np.any(idx3):
        indices = np.floor(x_arr[idx3] / L).astype(int)
        for i, idx in enumerate(np.where(idx3)[0]):
            V[idx] = coulomb_vec(x_arr[idx], indices[i]*L)
            V[idx] += coulomb_vec(x_arr[idx], (indices[i]+1)*L)
    
    # Return scalar if input was scalar
    return V[0] if np.isscalar(x) else V



