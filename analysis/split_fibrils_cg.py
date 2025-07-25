'''split aSyn fibril to independent protein chains'''
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import sys

traj = md.load_dcd('traj.dcd',top='top.pdb')

n_chains_aSyn = 209
N_res_aSyn = 140

top = md.Topology()
for _ in range(n_chains_aSyn):
    chain = top.add_chain()
    for residue in traj.top.chain(0).residues:
        if residue.index < N_res_aSyn:
            res = top.add_residue(residue.name, chain)      
            top.add_atom('CA', element=md.element.carbon, residue=res)

for i in range(1,traj.n_chains):
    chain = top.add_chain()
    for res in traj.top.chain(i).residues:
        res = top.add_residue(res.name, chain)
        top.add_atom('CA', element=md.element.carbon, residue=res)
        
new_traj = md.Trajectory(traj.xyz, top, traj.time, traj.unitcell_lengths, traj.unitcell_angles)

new_traj.save('new_traj.dcd')
new_traj[0].save('new_top.pdb')