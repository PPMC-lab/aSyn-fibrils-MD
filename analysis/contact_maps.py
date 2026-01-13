import numpy as np
import mdtraj as md
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--start',nargs='?',required=True,type=int)
parser.add_argument('--end',nargs='?',required=True,type=int)
args = parser.parse_args()
start=args.start
start_str=str(start)
end=args.end

# Functions to calculate energies
# Ashbaugh-Hatch potential
HALR = lambda r,s,l : 4*0.8368*l*((s/r)**12-(s/r)**6)
HASR = lambda r,s,l : 4*0.8368*((s/r)**12-(s/r)**6)+0.8368*(1-l)
HA = lambda r,s,l : np.where(r<2**(1/6)*s, HASR(r,s,l), HALR(r,s,l))
HASP = lambda r,s,l,rc : np.where(r<rc, HA(r,s,l)-HA(rc,s,l), 0)

# Debye-Hückel potential
DH = lambda r,yukawa_eps,lD : yukawa_eps*np.exp(-r/lD)/r
DHSP = lambda r,yukawa_eps,lD,rc : np.where(r<rc, DH(r,yukawa_eps,lD)-DH(rc,yukawa_eps,lD), 0)

def calc_contact_map(name,replica,N_1,N_2,temp,ionic,pH, start_str, start, end):
    """Calculate contact and energy maps between two proteins of different type
    in a system with equal number of chains of type A and B"""
    # load traj
    t = md.load(f'new_traj_{replica:s}.dcd',
        top=f'new_top_{replica:s}.pdb')[start:end] #Split traj into chunks, discard the first 100ns


    print(name,replica,t.n_frames, start, end)

    # load per-residue parameters
    df_residues = pd.read_csv(f'../input/residues_CALVADOS3.csv',index_col=1)

    # calculate traj of chain COM
    cmtop = md.Topology()
    xyz = np.empty((t.n_frames,t.n_chains,3))
    for chain in t.top.chains:
        atom_name = 'A' if chain.n_atoms == N_1 else 'B'
        #atom_name = 'A' if chain.n_atoms != N_2 else 'B' #all peptides
        mws = np.asarray([df_residues.loc[res.name].MW for res in chain.residues])

        new_chain = cmtop.add_chain()
        res = cmtop.add_residue('COM', new_chain, resSeq=chain.index)
        cmtop.add_atom(atom_name, element=t.top.atom(0).element, residue=res)
        t_chain = t.atom_slice(t.top.select(f'chainid {chain.index:d}'))
        com = np.sum(t_chain.xyz*mws[np.newaxis,:,np.newaxis],axis=1)/mws.sum()
        xyz[:,chain.index] = com
    cmtraj = md.Trajectory(xyz, cmtop, t.time, t.unitcell_lengths, t.unitcell_angles)

    L = t.unitcell_lengths[0][1]/2

    chains_1 = cmtraj.top.select('name A')
    chains_2 = cmtraj.top.select('name B')

    pairs = cmtraj.top.select_pairs('name A','name B')

    # calculate COM-COM separations and plot it
    dist = md.compute_distances(cmtraj,pairs)

    ###for minimum distances
    #dist = dist.reshape((-1,chains_1.size,chains_2.size))
    #plt.plot(np.min( np.min(dist, axis=2), axis=1 ))
    #plt.xlabel("# frame")
    #plt.ylabel("min COM-COM distance (nm)")
    #plt.title("Peptide-aSyn COM-COM")
    #plt.savefig('com_LL37_8a9l_min.png', dpi=300)
    #sys.exit()

    # calculate His charge
    df_residues.loc['H','q'] = 1. / ( 1 + 10**(pH-6) )

    qs = []
    lambdas = []
    sigmas = []
    fastas = []

    for chain in [t.top.chain(0),t.top.chain(-1)]:
        fasta = [res.name for res in chain.residues]
        fastas.append(fasta)
        sigmas.append(df_residues.loc[fasta].sigmas)
        lambdas.append(df_residues.loc[fasta].lambdas)
        charges = df_residues.loc[fasta].q.values
        charges[0] += 1
        charges[-1] -= 1
        qs.append(charges)

    sigma_pairs = 0.5*np.sum(np.asarray(list(itertools.product(sigmas[0],sigmas[1]))),axis=1)
    lambda_pairs = 0.5*np.sum(np.asarray(list(itertools.product(lambdas[0],lambdas[1]))),axis=1)
    q_pairs = np.prod(np.asarray(list(itertools.product(qs[0],qs[1]))),axis=1)

    RT = 8.3145*temp*1e-3
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/RT
    yukawa_eps_pairs = q_pairs*lB*RT
    lD = 1. / np.sqrt(8*np.pi*lB*ionic*6.022/10)

    #defining arrays
    tot_contact_pep = np.zeros((t.n_frames,chains_1.size)) #total contacts of each peptide with any asyn
    tot_contact_asyn = np.zeros(chains_2.size) #total contacts of each asyn with any peptide
    tot_contact_pep_asyn = np.zeros((chains_1.size,chains_2.size)) #total contacts of each peptide with each asyn chain
    dmap_contact = np.zeros((N_1,N_2)) #total residue contacts of each peptide and each asyn
    #dmap_contact_min = np.empty((chains_1.size*chains_2.size,N_1,N_2)) #initialize empty list to collect matrices
    #dmap_contact_min = np.full((chains_1.size*chains_2.size,N_1,N_2),  np.nan) #initialize empty list with nan to collect matrices
    
    for pair_ndx,(i,j) in enumerate(itertools.product(chains_1,chains_2)):
        #axis=0 is time, axis=1 is peptides, axis=2 is aSyn
        print(i,j,flush=True)
        frames_sel = dist[:,pair_ndx]<15 #calculate resiude distances only for selected COM pairs
        pair_indices = t.top.select_pairs(f'chainid {i:d}',f'chainid {j:d}')
        d = md.compute_distances(t[frames_sel],pair_indices)

        tot_contact_pep[frames_sel,i-chains_2.size] += (.5-.5*np.tanh((d-1.)/.3)).sum(axis=1) #continuos function to calculate contacts
        #tot_contact_asyn[j] +=  (.5-.5*np.tanh((d-1.)/.3)).sum()
        #tot_contact_pep_asyn[i-chains_2.size, j] += (.5-.5*np.tanh((d-1.)/.3)).sum()
        #dmap_contact += (.5-.5*np.tanh((d-1.)/.3)).sum(axis=0).reshape((N_1,N_2))
        
        #if d.size > 0:
        #    dmap_contact_min[pair_ndx] = d.min(axis=0).reshape((N_1,N_2)) #substitute each pair with minimum distance

    np.save(f'data_{replica:s}/{replica:s}_tot_contact_pep_chunk_{start_str:s}.npy',tot_contact_pep) #total contacts of each peptide with any asyn
    #np.save(f'data_{replica:s}/{replica:s}_tot_contact_asyn_chunk_{start_str:s}.npy',tot_contact_asyn) #avg contacts of each asyn with any peptide
    #np.save(f'data_{replica:s}/{replica:s}_tot_contact_pep_asyn_chunk_{start_str:s}.npy',tot_contact_pep_asyn) #avg contacts of each peptide with each asyn chain
    #np.save(f'data_{replica:s}/{replica:s}_contact_map_chunk_{start_str:s}.npy',dmap_contact) #avg residue contacts of each peptide with fiber
    #dmap_contact_min_absolute = np.nanmin(dmap_contact_min,axis=0)
    #dmap_contact_min_mean = np.nanmean(dmap_contact_min,axis=0)


# solution conditions
temp = 293
ionic = 0.1
pH = 7.5

N_res_dict = {'aSyn': 140, 'LL37': 37} #residues asyn and peptide

replica = '8a9l_180_LL37' #system_name

if not os.path.isdir(f'data_{replica:s}'):
    os.mkdir(f'data_{replica:s}')

calc_contact_map('LL37',replica,N_res_dict['LL37'],N_res_dict['aSyn'],temp,ionic,pH, start_str, start, end)
