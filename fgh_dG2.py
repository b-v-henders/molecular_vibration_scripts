"""
  Fourier Grid Hamiltonian method for 1-dimensional bound-state eigenvalues/eigenvectors.
  
  Original code:
  Tuguldur T. Odbadrakh
  Department of Chemistry
  University of Pittsburgh

  Modifications:
  Bryan Henderson
  Department of Chemistry
  University of Pittsburgh

  Reference: 
    D. J. Tannor, Introduction to Quantum Mechanics: A Time-Dependent Perspective, 2007, pp.297
    J. Stare, G. G. Balint-Kurti, J. Phys. Chem. A 2003, 107, 7204-7214
    M. A. Harthcock, J. Laane, J. Mol. Spec. 1982, 91, 300-324

  Usage:
    Import fgh_utils into the python script and call the function solve_FGH() with the following
    inputs:
      
      L - the length of the x axis
      V - array of the potential energy evaluated along the x axis
      G - Value of the element g44 of the G matrix (corresponding to 1/reduced mass) as a function of vibrational coordinate x
      dG2 - Second derivitive of G with respect to x

    The function will return the eigenvalues and eigenvectors of the resulting Hamiltonian in
    atomic units, i.e. hbar = 1.

"""

import numpy as np
from math import sin, cos, tan

def construct_H(x,V,G,dG2):
  """ Constructs the Hamiltonian matrix in atomic units. """

  N = np.size(V)
  H = np.zeros((N,N))
  pi = np.pi
  K = pi/((x[-1]-x[0])/N)
  Kb = 2.*pi/(x[-1]-x[0])
  n = (N-1.)/2.

  for i in range(N):
    for j in range(i+1):
      if i == j:
        H[i,j] = (G[i]/(2.))*(K**2.)*(1.0-1/(N**2.))/3. +0.25*dG2[i]+V[i]
      else:
        arg = pi*(j-i)/N
        H[i,j] = ((G[i] + G[j])/(4.))*(2.*(K**2.)*((-1)**(j-i))*cos(arg))/((N**2.)*sin(arg)*sin(arg))
        H[j,i] = H[i,j]
  return H

def solve_FGH(x,V,G,dG2):
  """ Solves the S.E. with potential V using the FGH method. """
  H = construct_H(x,V,G,dG2)
  return np.linalg.eigh(H)
