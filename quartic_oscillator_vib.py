"""
Script to solve two vibrational potentials with the same reduced mass.
Designed for solving rock-stretch coupling in water-anion clusters such as
HCO2-(H2O) or NO3-(H2O). See the following references for a more thourough
description of the problem:

Myshakin, E. M.; Jordan, K. D.; Sibert, E. L.; Johnson, M. A. J. Chem. Phys. 2003, 119 (19), 10138–10145.
Heine, N.; Kratz, E. G.; Bergmann, R.; Schofield, D. P.; Asmis, K. R.; Jordan, K. D.; McCoy, A. B. J. Phys. Chem. A 2014, 118 (37), 8188–8197.

Need the file fgh_dG2.py to solve to solve vibrational problem with FGH-DVR.

Usage:
python3 quartic_oscillator_vib.py input_file.txt

Input file should have the following format:

title line
x	V1	V2	G

Where x is the vibrational coordinate, and V1, V2, and G are the values of the two potentials and the associated
G-matrix element at x. In this implementation, V1 and G are assumed to be symmetric. V1 and V2 should be in
wavenumbers, and G should be in atomic units (NOT AMU!). Units of x are irrelevant as long as it is consistent with
the calculation of G. X should be increasing.

G is a value associated with the reduced mass of a vibration. See the following reference for details:

Rush, D. J.; Wiberg, K. B. J. Phys. Chem. A 1997, 101 (17), 3143–3151.

This script will print the first several energy levels of each potential, followed by the overlap, relative intensity
(overlap squared), and change in energy between the first state of V1 and the first several states in V2.

3 files will be output. energy.out simply outputs the first several energy levels for both potentials in cm-1 and
hartrees. transitions.out gives aformentioned transition energies and intentsities between the ground state of V1
and various states of V2. wavefunctions.out gives the wavefunctions of a few specific states, useful for plotting.

Written by Bryan Henderson, University of Pittsburgh

"""

import fgh_dG2 as fgh
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as uspline
import sys

N = 1001


def ho_dvr(x,V,G,dG2):        
        # solving the eigenvalue equation
        E, psi = fgh.solve_FGH(x,V,G,dG2)
        # print some stuff
        print("\nN = %i\n"%N)
        for i in range(17):
            print("%3i%15.8f%15.8f"%(i,E[i],np.sum(psi[:,i]**2.)))

        return E, psi

# HO parameters in a.u.
x = np.linspace(-60,60,N)

with open(sys.argv[1]) as f:
    pot_mass_file = f.readlines()
    f.close()

pot_mass = np.zeros((4,len(pot_mass_file) - 1))
for i in range (1,len(pot_mass_file)):
    line = pot_mass_file[i].split()
    for j in range (4):
        pot_mass[j,i-1]=float(line[j])


V1_polyfit = np.polynomial.polynomial.polyfit(pot_mass[0,:],pot_mass[1,:],np.array([0,2,4]))
V1= (V1_polyfit[2] * x**2 + V1_polyfit[4] * x**4) * 4.556335 * 10**-6

V2_polyfit = np.polynomial.polynomial.polyfit(pot_mass[0,:],pot_mass[2,:],np.array([0,1,2,3,4]))
V2 = (V2_polyfit[1] * x + V2_polyfit[2] * x**2 + V2_polyfit[3] * x**3 + V2_polyfit[4] * x**4 +V2_polyfit[0] - V1_polyfit[0]) * 4.556335 * 10**-6


u_polyfit = np.polynomial.polynomial.polyfit(pot_mass[0,:],1/pot_mass[3,:],np.array([0,2]))
u=u_polyfit[0]+u_polyfit[2]*x**2#+u_polyfit[4]*x**4

G1=1/u
du=2*u_polyfit[2]*x#+4*u_polyfit[4]*x**3
du2=2*u_polyfit[2]#+12*u_polyfit[4]*x**2

dG2=np.zeros((N))
for i in range(N):
    dG2[i] = 2*du[i]**2/(u[i]**3)-du2/(u[i]**2)


E1, psi1 = ho_dvr(x,V1,G1,dG2)
E2, psi2 = ho_dvr(x,V2,G1,dG2)


f = open("energy.out","w")
f.write("state n=0(hartree) n=1(hartree) n=0(cm-1) n=1(cm-1)\n")
for i in range(16):
    f.write("{0} {1} {2} {3} {4}\n".format(i, E1[i], E2[i], E1[i]*219474.63, E2[i]*219474.63))
f.close()
print()
print("Overlap of eigenstates:")
print()
print("state    overlap        intensity		Delta E")

f = open("transition.out","w")
f.write("state deltaE(cm-1) overlap intensity\n")
for itera in range (0,16):
    overlap = 0
    intensity = 0
    exp_x = 0
    for iterb in range (0,N):   
        overlap += psi1[iterb,0] * psi2[iterb,itera]
        #exp_x += psi2[iterb,itera]**2 * x[iterb]
    intensity = overlap**2.
    deltaE = (E2[itera]-E1[0]) * 219474.63
    print(str(itera) + "     " + str(overlap) + "      " + str(intensity) + "       " + str(deltaE))
    f.write("{0} {1} {2} {3}\n".format(itera, deltaE, overlap, intensity))
f.close()


n0_levels=1
n1_levels=4
with open("wavefunctions.out","w") as f:
	f.write("{0}   {1}\n".format(n0_levels,n1_levels))
	f.write("x")
	for i in range(n0_levels):
		f.write("   N=0_n={0}".format(i))
	for i in range(n1_levels):
		f.write("   N=1_n={0}".format(i))
	f.write("\n")
	
	f.write("null")
	for i in range(n0_levels):
		f.write("   {0}".format(E1[i]/(4.556335 * 10**-6)))
	for i in range(n1_levels):
		f.write("   {0}".format(E2[i]/(4.556335 * 10**-6)))
	f.write("\n")
	for i in range(len(x)):
		f.write(str(x[i]))
		
		for j in range(n0_levels):
			f.write("    {0}".format(psi1[i,j]))
		for j in range(n1_levels):
			f.write("    {0}".format(psi2[i,j]))
		f.write("\n")
