import numpy as np
import mdtraj as md
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import os
import sys
import scipy.stats as scs
sys.path.append('BLOCKING')
from main import BlockAnalysis

def autoblock(cv, multi=1):
	'''autoblock function to calculate the blocking error'''
	block = BlockAnalysis(cv, multi=multi)
	block.SEM()
	return block.sem, block.bs

def kde(a):
    min_ = np.min(a)
    max_ = np.max(a)
    x = np.linspace( min_, max_, num = 100 )
    d = scs.gaussian_kde( a, bw_method = "silverman" ).evaluate(x)
    u = np.average(a)
    return x,d/np.sum(d),u

#traj = md.load_dcd('new_traj_pacap_cubic.dcd','new_top_pacap_cubic.pdb')[1000:]
traj = md.load_dcd('new_traj_pacap_alone.dcd','new_top_pacap_alone.pdb')[1000:]

df_residues = pd.read_csv(f'../input/residues_CALVADOS3.csv',index_col=0)
fasta="LGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA" #C-t
fasta_pep="HSDGIFTDSYSRYRKQMAVKKYLAAVLGKRYKQRVKNK" #pacap

masses = df_residues.loc[list(fasta_pep),'MW'].values
masses[0] += 2
masses[-1] += 16

n_chains=105
n=1
i=0 #start C-t
j=37 #end C-t
 
rg_asyn_chains=[]
hist_rg_asyn_chains=[]
hist_rg_prob_asyn_chains=[]

for _ in range(n_chains): #calculate Rg for all aSyn chains in fiber
	print("Chain_"+str(n))
	print(i,j)
	traj_idx = traj.atom_slice(traj.top.select(f'index {i:d} to {j:d}'))


	# calculate the center of mass
	cm = np.sum(traj_idx.xyz*masses[np.newaxis,:,np.newaxis],axis=1)/masses.sum()

	# calculate residue-cm distances
	si = np.linalg.norm(traj_idx.xyz - cm[:,np.newaxis,:],axis=2)

	# calculate rg
	rg_array = np.sqrt(np.sum(si**2*masses,axis=1)/masses.sum())
	rg_hist=kde(rg_array)


	rg_mean = np.mean(rg_array)
	rg_se, rg_blocksize = autoblock(rg_array)
	
	rg_asyn_chains.append(rg_mean)
	#hist_rg_asyn_chains.append(rg_hist[0])
	#hist_rg_prob_asyn_chains.append(rg_hist[1])

	i+=38
	j+=38
	n+=1
	
	#plt.plot(rg_hist[0],rg_hist[1], 'b')
	#plt.plot(rg_hist_alone[0],rg_hist_alone[1], 'r')
	#plt.show()



print(rg_asyn_chains)

np.save(f'./data_rg/rg_pacap_alone.npy',rg_asyn_chains)