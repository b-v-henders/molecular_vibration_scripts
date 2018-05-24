"""

Let's combine all of the tedious geometry stuff I needed to do into one program.

Usage:
python3 G_matrix_all.py input

Input file format:

num_atoms	num_geometries
output_filename
atom1_name	mass_atom1
atom2_name	mass_atom2
...
atom_n_name	mass_atom_n

geom1_name
atom1_x	atom1_y	atom1_z
...
atom_n_x	atom_n_y	atom_n_z

geom2_name
atom1_x	atom1_y	atom1_z
...


All geometries should be in angstroms, and all masses should be in atomic units (me) (NOT AMU!!!).
Translates coordinates into center of mass coordinates in bohrs, rotates within the Eckart frame,
calculates derivatives, and calculates the vibrational g matrix term for each geometry.

First geometry is assumed to be the reference geometry for Eckart frame stuff.

Vibrational coordinate can be defined in the function coord_define. In this case, the coordinate
is the rock angle in a water-ion cluster. This is system specific, and requires modification for
different systems.

This script assumes symmetry around the vibrational coordinate. This is explicitly difined in the 
function reflect. Again, this is system specific.

The Eckart frame rotation only occurs in two dimensions (XY). This should be valid as long as there's
a plane of symmetry in the Z plane that is retained throughout the motion. If your molecule is planar,
ensure it is in the XY plane. Also, make sure your coordinates are consistent throughout the motion
(ie no unexpected changes in sign, swapping x and y coordinates, etc.).

Output file contains all transformed geometries (in bohrs), and derivatives, vibrational coordinate,
g-matrix element, and reduced mass for each point.

Written by Bryan Henderson, University of Pittsburgh

"""

import sys
import numpy as np
import math
from scipy.interpolate import UnivariateSpline as uspline

#This function defines your vibrational coordinate. System specific, modify as necessary.
#Make sure this is unique for each geometry
def coord_define(geom):
	O_H1_dist=atom_dist(geom[0,:],geom[1,:])
	O_H2_dist=atom_dist(geom[0,:],geom[2,:])
	H_H_dist=atom_dist(geom[1,:],geom[2,:])
	HOH_angle=law_cosines(O_H1_dist,O_H2_dist,H_H_dist)
	
	O_C_dist=atom_dist(geom[0,:],geom[3,:])
	H1_C_dist=atom_dist(geom[1,:],geom[3,:])
	H1_O_C_angle=law_cosines(O_H1_dist,O_C_dist,H1_C_dist)
	
	phi=H1_O_C_angle-0.5*HOH_angle
	if abs(phi)<10**-7:
		return 0.
	else:
		return phi

#reflects the geometry for positive and negative coordinates. System specific, modify as necessary.
def reflect(geom_ref):
	geom_reflect=np.zeros((natoms,3))
	for i in range(natoms):
		if i==1:
			geom_reflect[2,0]=geom_ref[1,0]
			geom_reflect[2,1]=-1*geom_ref[1,1]
			geom_reflect[2,2]=geom_ref[1,2]
		elif i==2:
			geom_reflect[1,0]=geom_ref[2,0]
			geom_reflect[1,1]=-1*geom_ref[2,1]
			geom_reflect[1,2]=geom_ref[2,2]
		elif i==4:
			geom_reflect[6,0]=geom_ref[4,0]
			geom_reflect[6,1]=-1*geom_ref[4,1]
			geom_reflect[6,2]=geom_ref[4,2]
		elif i==6:
			geom_reflect[4,0]=geom_ref[6,0]
			geom_reflect[4,1]=-1*geom_ref[6,1]
			geom_reflect[4,2]=geom_ref[6,2]
		else:
			geom_reflect[i,0]=geom_ref[i,0]
			geom_reflect[i,1]=-1*geom_ref[i,1]
			geom_reflect[i,2]=geom_ref[i,2]
	return geom_reflect

def derv_spl(coordinate,position):
	spl=uspline(coordinate,position,k=5)
	dervs=np.zeros((len(coordinate)))
	for i in range(len(coordinate)):
		dervs[i]=spl.derivatives(coordinate[i])[1]
	return(dervs)

def atom_dist(atom1_pos,atom2_pos):
	distance_sq=0.
	for i in range(3):
		distance_sq+=(atom1_pos[i]-atom2_pos[i])**2
	return np.sqrt(distance_sq)

def law_cosines(dist1,dist2,dist3):
	return (np.arccos((dist1**2+dist2**2-dist3**2)/(2*dist1*dist2))*(180/np.pi))

def eck_value(ref,geom,mass,atoms):
    eck_rot_value = 0
    for i in range(atoms):
        eck_rot_value += mass[i]*(ref[i,0] * (geom[i,1]-ref[i,1]) - ref[i,1] * (geom[i,0]-ref[i,0]))
    return eck_rot_value

def rotate(geom,theta,atoms):
    new_geom = np.zeros((atoms,2))
    for i in range(atoms):
        new_geom[i,0] = math.cos(theta)*geom[i,0] + math.sin(theta)*geom[i,1]
        new_geom[i,1] = -1 * math.sin(theta) * geom[i,0] + math.cos(theta) * geom[i,1]
    return new_geom

#This bit is an absolute mess that I don't know how it works.
#Making this not terrible should be done at some point.
def eckart_rotate(reference_geom,current_geom,mass,num_atoms):
    largestep = 0.00001
    smallstep = 0.00000001
    max_iterations_l = 100000
    max_iterations_s = int(largestep/smallstep)
    
    eck_rot_old = eck_value(reference_geom,current_geom,mass,num_atoms)
    test_geom = rotate(current_geom,largestep,num_atoms)
    eck_rot_new = eck_value(reference_geom,test_geom,mass,num_atoms)
    if abs(eck_rot_new)<abs(eck_rot_old):
        for j in range(2,max_iterations_l):
            eck_rot_old = eck_rot_new
            test_geom = rotate(current_geom,largestep * j,num_atoms)
            eck_rot_new = eck_value(reference_geom,test_geom,mass,num_atoms)
            if abs(eck_rot_new)>abs(eck_rot_old):
                eck_rot_old = eck_rot_new
                test_geom_old = test_geom
                test_geom=rotate(current_geom,largestep * (j-1) + smallstep,num_atoms)
                eck_rot_new = eck_value(reference_geom,test_geom,mass,num_atoms)
                if abs(eck_rot_new)<abs(eck_rot_old):
                    for k in range(2,max_iterations_s):
                        eck_rot_old = eck_rot_new
                        test_geom_old = test_geom
                        test_geom = rotate(reference_geom,largestep * (j-1) + smallstep * k,num_atoms)
                        eck_rot_new = eck_value(reference_geom,test_geom,mass,num_atoms)
                        if abs(eck_rot_new)>abs(eck_rot_old):
                            converge=1
                            break
                else:
                    for k in range(1,max_iterations_s):
                        eck_rot_old = eck_rot_new
                        test_geom_old = test_geom
                        test_geom = rotate(reference_geom,largestep * (j-1) - smallstep * k,num_atoms)
                        eck_rot_new = eck_value(reference_geom,test_geom,mass,num_atoms)
                        if abs(eck_rot_new)>abs(eck_rot_old):
                            converge=1
                            break
                break
            if j == max_iterations_l - 1:
                print("no convergance")
                converge=0
    else:
        for j in range(1,max_iterations_l):
            eck_rot_old = eck_rot_new
            test_geom = rotate(current_geom,-1*largestep * j,num_atoms)
            eck_rot_new = eck_value(reference_geom,test_geom,mass,num_atoms)
            if abs(eck_rot_new)>abs(eck_rot_old):
                eck_rot_old = eck_rot_new
                test_geom_old = test_geom
                test_geom=rotate(current_geom,-1*largestep * (j-1) + smallstep,num_atoms)
                eck_rot_new = eck_value(reference_geom,test_geom,mass,num_atoms)
                if abs(eck_rot_new)<abs(eck_rot_old):
                    for k in range(2,max_iterations_s):
                        eck_rot_old = eck_rot_new
                        test_geom_old = test_geom
                        test_geom = rotate(reference_geom,-1*largestep * (j-1) + smallstep * k,num_atoms)
                        eck_rot_new = eck_value(reference_geom,test_geom,mass,num_atoms)
                        if abs(eck_rot_new)>abs(eck_rot_old):
                            converge=1
                            break
                else:
                    for k in range(1,max_iterations_s):
                        eck_rot_old = eck_rot_new
                        test_geom_old = test_geom
                        test_geom = rotate(reference_geom,-1*largestep * (j-1) - smallstep * k,num_atoms)
                        eck_rot_new = eck_value(reference_geom,test_geom,mass,num_atoms)
                        if abs(eck_rot_new)>abs(eck_rot_old):
                            converge=1
                            break
                break
            if j == max_iterations_l - 1:
                print("no convergance")
                converge=0
    return test_geom_old

def find_gmat(geom,derivative,mass,natoms):
	gmat_inv=np.zeros((4,4))
	
	#pure rotational terms
	for coord1 in range(3):
		for coord2 in range(coord1,3):
			Iik=0
			if coord1 == coord2:
				for atom in range(natoms):
					Iik+=mass[atom]*(geom[atom,0]**2+geom[atom,1]**2+geom[atom,2]**2-geom[atom,coord1]**2)
			else:
				for atom in range(natoms):
					Iik+=-1*mass[atom]*geom[atom,coord1]*geom[atom,coord2]
			gmat_inv[coord1,coord2]=Iik
			gmat_inv[coord2,coord1]=Iik
	
	#rotational-vibrational terms
	for coord in range(3):
		Xi=0
		for atom in range(natoms):
			Xi+=mass[atom]*np.cross(geom[atom,:],derivative[atom,:])[coord]
		gmat_inv[coord,3]=Xi
		gmat_inv[3,coord]=Xi
	
	#pure vibraional term
	Yi=0
	for coord in range(3):
		for atom in range(natoms):
			Yi+=mass[atom]*derivative[atom,coord]**2
	gmat_inv[3,3]=Yi
	
	gmat=np.linalg.inv(gmat_inv)
	
	return gmat_inv,gmat
			

#Coordinate file in Angstroms
with open(sys.argv[1]) as f:
	coordinate_file=f.readlines()
	f.close()

natoms=int(coordinate_file[0].split()[0])
ngeom=int(coordinate_file[0].split()[1])

masses=np.zeros((natoms))
atomnames=[]

for i in range(natoms):
	line=coordinate_file[i+2].split()
	atomnames.append(line[0])
	masses[i]=line[1]

geom_original=np.zeros((ngeom,natoms,3))
geom_name=[]

for geom in range(ngeom):
	geom_name.append(coordinate_file[10+(2+natoms)*geom].split()[0])
	for atom in range(natoms):
		line=coordinate_file[11+(2+natoms)*geom+atom].split()
		for coord in range(3):
			geom_original[geom,atom,coord]=line[coord]
	
geom_com_bohr=np.zeros((ngeom,natoms,3))

mass_tot=np.sum(masses[:])
print("Converting to center of mass coordinates in bohrs...")
for geom in range(ngeom):
	for coord in range(3):
		com=0
		for atom in range(natoms):
			com+=geom_original[geom,atom,coord]*masses[atom]/mass_tot
		
		for atom in range(natoms):
			geom_com_bohr[geom,atom,coord]=(geom_original[geom,atom,coord]-com)/0.529177249

print("Conversion complete!")
print("Rotating system with respect to Eckart frame...")
geom_eck=np.zeros((ngeom,natoms,3))

geom_eck[0,:,:]=geom_com_bohr[0,:,:]
geom_eck[:,:,2]=geom_com_bohr[:,:,2]


for geom in range(1,ngeom):
	geom_eck[geom,:,0:2]=eckart_rotate(geom_eck[0,:,0:2],geom_com_bohr[geom,:,0:2],masses,natoms)
print("Rotation complete!")	
#print(geom_eck)

print("Defining vibrational coordinate system...")

vib_coord=np.zeros((ngeom))

for i in range(ngeom):
	vib_coord[i]=coord_define(geom_eck[i,:,:])

vib_order=np.argsort(vib_coord)

ngeom_fin = 2*ngeom-1
geom_fin=np.zeros((ngeom_fin,natoms,3))

vib_coord_fin=np.zeros((ngeom_fin))

for i in range(ngeom):
	if vib_coord[i]==0:
		vib_coord_fin[ngeom-1]=vib_coord[i]
		geom_fin[ngeom-1,:,:]=geom_eck[i,:,:]
	elif vib_coord[i]>0:
		vib_coord_fin[ngeom+vib_order[i]-1]=vib_coord[i]
		vib_coord_fin[ngeom-vib_order[i]-1]=-1*vib_coord[i]
		geom_fin[ngeom+vib_order[i]-1,:,:]=geom_eck[i,:,:]
		geom_fin[ngeom-vib_order[i]-1,:,:]=reflect(geom_eck[i,:,:])
	else:
		vib_coord_fin[ngeom+vib_order[i]]=-1*vib_coord[i]
		vib_coord_fin[ngeom-vib_order[i]-1]=vib_coord[i]
		geom_fin[ngeom+vib_order[i]-1,:,:]=reflect(geom_eck[i,:,:])
		geom_fin[ngeom-vib_order[i]-1,:,:]=geom_eck[i,:,:]

print("Coordinate system definition complete!")
print("Calculating derivatives...")

derivatives=np.zeros((ngeom_fin,natoms,3))

for atom in range(natoms):
	for coord in range(3):
		derivatives[:,atom,coord]=derv_spl(vib_coord_fin,geom_fin[:,atom,coord])

print("Derivatives complete!")
print("Calculating G-matrix...")
G_inverse=np.zeros((ngeom_fin,4,4))
G_mat=np.zeros((ngeom_fin,4,4))
for i in range((ngeom_fin)):
	G_inverse[i,:,:],G_mat[i,:,:]=find_gmat(geom_fin[i,:,:],derivatives[i,:,:],masses,natoms)

print("G-matrix calculated!")

print("Generating output...")
with open(coordinate_file[1].split()[0],'w') as f:
	f.write("Atom	Mass	X	Y	Z	dX/dq	dY/dq	dZ/dq	Phi	G	Reduced_mass\n")
	for geom in range(ngeom_fin):
		f.write("{0}	{1}	{2}	{3}	{4}	{5}	{6}	{7}	{8}	{9}	{10}\n".format(atomnames[0],masses[0],geom_fin[geom,0,0],geom_fin[geom,0,1],geom_fin[geom,0,2],derivatives[geom,0,0],derivatives[geom,0,1],derivatives[geom,0,2],vib_coord_fin[geom],G_mat[geom,3,3],1/G_mat[geom,3,3]))
		for atom in range(1,natoms):
			f.write("{0}	{1}	{2}	{3}	{4}	{5}	{6}	{7}\n".format(atomnames[atom],masses[atom],geom_fin[geom,atom,0],geom_fin[geom,atom,1],geom_fin[geom,atom,2],derivatives[geom,atom,0],derivatives[geom,atom,1],derivatives[geom,atom,2]))
		f.write("\n\n")

print("Everything complete!")













