import numpy as np
import sys

def index_ijk(a,b,c):
	if a>=b:
		if b>=c:
			i=a+1
			j=b+1
			k=c+1
		elif a>=c:
			i=a+1
			j=c+1
			k=b+1
		else:
			i=c+1
			j=a+1
			k=b+1
	else:
		if a>=c:
			i=b+1
			j=a+1
			k=c+1
		elif c>=b:
			i=c+1
			j=b+1
			k=a+1
		else:
			i=b+1
			j=c+1
			k=a+1
	
	index=(1/6)*(i**3-i)+0.5*(j**2-j)+k-1
	return int(index)
  
def index_ijkl(a,b,c,d):
	indicies=np.zeros((4))
	indicies[0]=a+1
	indicies[1]=b+1
	indicies[2]=c+1
	indicies[3]=d+1
	
	indicies_sorted = np.sort(indicies)
	i=indicies_sorted[3]
	j=indicies_sorted[2]
	k=indicies_sorted[1]
	l=indicies_sorted[0]
	index=(1/24)*(i**4+2*i**3-i**2-2*i)+(1/6)*(j**3-j)+0.5*(k**2-k)+l-1
	return int(index)
	

def quanta_greater(n,m):
	if n >= m:
		return n
	else:
		return m

def raising(ho_state,mode):
	ho_state_new=ho_state.copy()
	ho_state_new[mode,1]=np.sqrt(ho_state[mode,0]+1)*ho_state[mode,1]
	ho_state_new[mode,0] = ho_state[mode,0]+1
	return ho_state_new

def lowering(ho_state,mode):
	ho_state_new=ho_state.copy()
	ho_state_new[mode,1]=np.sqrt(ho_state[mode,0])*ho_state[mode,1]
	ho_state_new[mode,0]=ho_state[mode,0]-1
	if ho_state_new[mode,0]<0:
		ho_state_new[mode,0]=0
	return ho_state_new

def readquanta(ci_state,ci_index,nmodes):
	ho_state = np.ones((nmodes,2))
	for i in range(nmodes):
		ho_state[i,0]=ci_index[ci_state,i]
	return ho_state

def check_couple(couple_matrix,ho_state_1, ho_state_2):
	couple_tf=1
	for i in range(len(couple_matrix)):
		order_i=abs(ho_state_1[i,0]-ho_state_2[i,0])
		if couple_matrix[i]==0:
			if order_i != 0:
				couple_tf = 0
				break
		elif couple_matrix[i]==1:
			if order_i !=1:
				couple_tf = 0
				break				
		elif couple_matrix[i]==2:
			if order_i != 2 and order_i != 0:
				couple_tf = 0
				break
		elif couple_matrix[i]==3:
			if order_i != 3 and order_i != 1:
				couple_tf = 0
				break
		elif couple_matrix[i]==4:
			if order_i != 4 and order_i != 2 and order_i != 0:
				couple_tf = 0
				break
				
	return(couple_tf)

def integral(ho_states_i,ho_states_j,nmodes):
	integrand=1.
	for mode in range(nmodes):
		if ho_states_i[mode,0]!=ho_states_j[mode,0]:
			integrand = 0.
			break
		else:
			integrand*=ho_states_i[mode,1]*ho_states_j[mode,1]
	return integrand

def coupling(ho_states_i,ho_states_j,couple_matrix,couple_constant,nmodes):
	mode_coupling=0.
	temp_couple_matrix=couple_matrix.copy()
	temp_ho_states_j = [ho_states_j.copy()]
	couple_order=int(sum(couple_matrix))

	for iteration in range(couple_order):
		for mode in range(nmodes):
			if temp_couple_matrix[mode]>0:
				states_length=len(temp_ho_states_j)
				for i in range(states_length):
					new_state_down=temp_ho_states_j[i].copy()
					new_state_down=(lowering(new_state_down,mode))
					new_state_up=temp_ho_states_j[i].copy()
					new_state_up=raising(new_state_up,mode)
					temp_ho_states_j.append(new_state_down.copy())
					temp_ho_states_j[i]=new_state_up.copy()
				temp_couple_matrix[mode]+=-1				
	
	x_constants =(1/(np.sqrt(2)))**couple_order				
			
	for i in range(len(temp_ho_states_j)):
		mode_coupling+=integral(ho_states_i,temp_ho_states_j[i],nmodes)

	return (mode_coupling*x_constants*couple_constant)

def momentum_squared(ho_states_i,ho_states_j,mode,frequency,nmodes):
	mod_j_up_up=raising(raising(ho_states_j.copy(),mode),mode)
	mod_j_up_down=lowering(raising(ho_states_j.copy(),mode),mode)
	mod_j_down_up=raising(lowering(ho_states_j.copy(),mode),mode)
	mod_j_down_down=lowering(lowering(ho_states_j.copy(),mode),mode)
	
	p_2_dimensionless=-1*integral(ho_states_i,mod_j_up_up,nmodes)-integral(ho_states_i,mod_j_down_down,nmodes)+integral(ho_states_i,mod_j_up_down,nmodes)+integral(ho_states_i,mod_j_down_up,nmodes)
	p_2=p_2_dimensionless*0.25*frequency
	return p_2

# CI input file
with open(sys.argv[1]) as f:
	ci_file=f.readlines()
	f.close()

nmodes = int(ci_file[0].split()[0])
ncubic = int(ci_file[0].split()[1])
ncubic_tot = index_ijk(nmodes-1,nmodes-1,nmodes-1)+1
nquartic = int(ci_file[0].split()[2])
nquartic_total=index_ijkl(nmodes-1,nmodes-1,nmodes-1,nmodes-1)+1

ci_index_filename=str(ci_file[1].split()[0])
ci_matrixfile_filename=str(ci_file[2].split()[0])
ci_output_filename=str(ci_file[3].split()[0])

freq = np.zeros((nmodes))
cubic = np.zeros((ncubic_tot))
quartic=np.zeros((nquartic_total))

for i in range(nmodes):
	freq[i]=float(ci_file[i+5])

for iteration in range(ncubic):
	line=ci_file[iteration+6+nmodes].split()
	i=int(line[0])
	j=int(line[1])
	k=int(line[2])
	cubic[index_ijk(i-1,j-1,k-1)]=float(line[3])

for iteration in range(nquartic):
	line=ci_file[iteration+7+nmodes+ncubic].split()
	i=int(line[0])-1
	j=int(line[1])-1
	k=int(line[2])-1
	l=int(line[3])-1
	quartic[index_ijkl(i,j,k,l)]=float(line[4])
			
# CI index file
with open(ci_index_filename)as f:
	ci_index_file = f.readlines()
	f.close

num_ci=len(ci_index_file)
ci_index=np.zeros((num_ci,nmodes+1))
for i in range(len(ci_index_file)):
	line=ci_index_file[i].split()
	order=0
	for j in range(nmodes):
		ci_index[i,j]=int(line[j])
		order+=int(line[j])
	ci_index[i,-1]=order

ci_matrix=np.zeros((num_ci,num_ci))

#HO states matrix: [i,0] is quanta in mode i, [i,1] is constant resulting from raising/lowering operators relative to mode i

for i in range(num_ci):
	states_i=readquanta(i,ci_index,nmodes)
	for j in range(num_ci):
		ci_element=0
		states_j=readquanta(j,ci_index,nmodes)

		index_difference = 0
		for k in range(nmodes):
			index_difference+=abs(states_i[k,0]-states_j[k,0])
		
		if index_difference < 5:
			if index_difference == 0 or index_difference == 2:
				for mode in range(nmodes):
					couple_matrix=np.zeros((nmodes))
					couple_matrix[mode]=2
					if check_couple(couple_matrix,states_i,states_j)==1:
					#momentum term
						ci_element+=momentum_squared(states_i,states_j,mode,freq[mode],nmodes)
					#second_order coupling
						coupling_matrix = np.zeros((nmodes))
						coupling_matrix[mode]=2
						ci_element+=coupling(states_i,states_j,coupling_matrix,freq[mode]/2,nmodes)
			
			if index_difference == 1 or index_difference == 3:
				#Third order coupling
				for mode1 in range(nmodes):
					for mode2 in range(nmodes):
						for mode3 in range(nmodes):
							coupling_matrix=np.zeros((nmodes))
							coupling_matrix[mode1]+=1
							coupling_matrix[mode2]+=1
							coupling_matrix[mode3]+=1
							if cubic[index_ijk(mode1,mode2,mode3)] != 0 and check_couple(coupling_matrix,states_i,states_j)==1:

								ci_element_cubic=coupling(states_i,states_j,coupling_matrix,cubic[index_ijk(mode1,mode2,mode3)]/6,nmodes)
								ci_element+=ci_element_cubic
								#print(str(mode1+1)+"   "+str(mode2+1)+"    "+str(mode3+1)+"   "+str(ci_element_cubic)+"   "+str(i)+"  "+str(j)+"   "+str(states_i)+"   "+str(states_j))
			if index_difference == 0 or index_difference == 2 or index_difference ==4:
				#Fourth order coupling
				for mode1 in range(nmodes):
					for mode2 in range(mode1,nmodes):
						for mode3 in range(mode2,nmodes):
							for mode4 in range(mode3,nmodes):
								coupling_matrix=np.zeros((nmodes))
								coupling_matrix[mode1]+=1
								coupling_matrix[mode2]+=1
								coupling_matrix[mode3]+=1
								coupling_matrix[mode4]+=1
								if quartic[index_ijkl(mode1,mode2,mode3,mode4)] != 0 and check_couple(coupling_matrix,states_i,states_j):

									ci_element_quartic=coupling(states_i,states_j,coupling_matrix,quartic[index_ijkl(mode1,mode2,mode3,mode4)]/24,nmodes)
									ci_element+=ci_element_quartic



		ci_matrix[i,j]=ci_element
				
	

#print(ci_matrix[1,:])

for i in range(num_ci):
	for j in range(num_ci):
		if abs(ci_matrix[i,j])<10**-10:
			ci_matrix[i,j]=0

E,psi=np.linalg.eigh(ci_matrix)

print(E)
E_excite=np.zeros((num_ci))
for i in range(num_ci):
	E_excite[i]=E[i]-E[0]
print(E_excite)
#print(psi[:,0])

#print(ci_index)
#print(freq)
#print(quartic_iiii)
#print(quartic_iijj)
with open(ci_matrixfile_filename,"w") as f:
	f.write("index    ")
	for i in range(num_ci):
		for j in range(nmodes):
			f.write("{0},".format(str(int(ci_index[i,j]))))
		f.write("    ")
	f.write("\n")
	
	for i in range(num_ci):
		for mode in range(nmodes):
			f.write("{0},".format(str(int(ci_index[i,mode]))))
		f.write("   ")
		
		for j in range(num_ci):
			f.write("{0}    ".format(ci_matrix[i,j]))
			#if ci_matrix[i,j] != ci_matrix[j,i]:
			#	print("Somethin's goofy    " + str(i) + "    " + str(j))
		f.write("\n")
	f.close()

with open("ci_wavefunction.txt","w") as f:
	f.write("index")
	for i in range(num_ci):
		f.write("    {0}".format(i))
	f.write("\n")
	
	for i in range(num_ci):
		for mode in range(nmodes):
			f.write("{0},".format(str(int(ci_index[i,mode]))))
		f.write("   ")
		
		for j in range(num_ci):
			f.write("{0}    ".format(psi[i,j]))
		f.write("\n")
	f.close()

with open(ci_output_filename,"w") as f:
	f.write("State  Total_Energy   Excitation_Energy\n")
	for i in range(num_ci):
		f.write("{0}    {1}    {2}".format(i,E[i],E_excite[i]))
		f.write("\n")
		
print(psi[:,0])
