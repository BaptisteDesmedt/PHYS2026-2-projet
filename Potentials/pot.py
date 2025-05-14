import numpy as np 
import math
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
    Q, e, a, K = params 
    
    if x < L:
        d = np.sqrt((x-L)**2 + a**2)
        V = (-(Q*e*K))/d
        dr = np.sqrt((x-(2*L))**2 + a**2)
        V += (-(Q*e*K))/dr
        return V
    if x > N*L:
        d = np.sqrt((x-N*L)**2 + a**2)
        V = (-(Q*e*K))/d
        dl = np.sqrt((x-((N-1)*L))**2 + a**2)
        V += (-(Q*e*K))/dl
        return V
    else : 
        indice = math.floor(x/L)
        dl = np.sqrt((x-(indice*L))**2 + a**2)
        V = -(Q*e*K)/(dl)
        dr = np.sqrt((x-((indice + 1)*L))**2 + a**2)
        V += -(Q*e*K)/(dr)
        return V
    return 0

def pot_solve(x,params,N,L):
    """
    Calculate the potential at a point or an array of points x due to a charge distribution.
    
    Parameters:
    x : float or ndarray
        Point(s) where the potential is calculated.
    params : tuple
        A tuple containing the parameters (Q, e, a, K)
    N : int
        Number of periods
    L : float
        Period length
    
    Returns:
    float or ndarray
        The potential at x.
    """
     
    Q, e, a, K = params  
    x_arr = np.atleast_1d(x)
    
    V = np.empty_like(x_arr, dtype=float)
    
    # Region where x < L
    idx1 = x_arr < L
    if np.any(idx1):
        d = ((x_arr[idx1] - L)**2 + a**2)**0.5
        V[idx1] = -(Q * e*K) / (d)
        dr = ((x_arr[idx1] - (2*L))**2 + a**2)**0.5
        V[idx1] += -(Q * e*K) / (dr)
    
    # Region where x > N*L
    idx2 = x_arr > N * L
    if np.any(idx2):
        d = ((x_arr[idx2] - N * L)**2 + a**2)**0.5
        V[idx2] = -(Q * e*K)/(d)
        dl = ((x_arr[idx2] - ((N-1)*L))**2 + a**2)**0.5
        V[idx2] += -(Q * e*K)/(dl)
    
    # Region in between
    idx3 = ~(idx1 | idx2)
    if np.any(idx3):
        # Use vectorized floor for indices
        indices = np.floor(x_arr[idx3] / L)
        dl = ((x_arr[idx3] - indices * L)**2 + a**2)**2
        dr = ((x_arr[idx3] - ((indices + 1) * L))**2 + a**2)**2
        V[idx3] = (-(Q * e*K)/(dl)) + (-(Q * e*K)/( dr))
    
    # Return scalar if input was scalar
    return V[0] if np.isscalar(x) or x_arr.shape == () else V



