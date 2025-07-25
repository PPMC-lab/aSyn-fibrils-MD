'''1 fibril as single MDP, multiple disordered peptides'''
import os
import pandas as pd
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--fiber',nargs='?',required=False,type=str)
parser.add_argument('--peptide',nargs='?',required=False,type=str)

args = parser.parse_args()

cwd = os.getcwd()

# set the side length of the cubic box and cylinder radius
L = 120

# set the saving interval (number of integration steps)
N_save = 10000

# set final number of frames to save
N_frames = 10000

#Simulation time (ns)
sim_time=str(int(N_save*N_frames*0.01/1000))

sysname = f'{args.fiber:s}_{args.peptide:s}_{sim_time:s}ns'
residues_file = f'{cwd}/input/residues_CALVADOS3.csv'
fasta_file = f'{cwd}/input/fastalib.fasta'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [L, L, L], # nm
  temp = 293.15, # K
  ionic = 0.15, # molar
  pH = 7.5,
  topol = 'random',

  # RUNTIME SETTINGS
  wfreq = N_save, # dcd writing interval, 1 = 10 fs
  steps = N_frames*N_save, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = 'CUDA', # or CUDA
  restart = 'checkpoint',
  frestart = 'restart.chk',
  verbose = True,

  # JOB SETTINGS (ignore if running locally)
  submit = False
)

# PATH
path = f'{cwd}/{sysname:s}'
subprocess.run(f'mkdir -p {path}',shell=True)
subprocess.run(f'mkdir -p data',shell=True)

analyses = f"""

from calvados.analysis import save_rg

save_rg("{path:s}","{sysname:s}","{residues_file:s}","data",10)
"""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = 1, # number of molecules
  restraint = True, # apply restraints
  ext_restraint = False, # apply external restraints
  charge_termini = 'both', # charge N or C or both
  
  # INPUT
  fresidues = residues_file, # residue definitions
  fdomains = f'{cwd}/input/domains.yaml', # domain definitions (harmonic restraints)
  pdb_folder = f'{cwd}/input', # directory for pdb and PAE files
  ffasta = fasta_file, # domain definitions (harmonic restraints)
  
  # RESTRAINTS
  restraint_type = 'harmonic', # harmonic or go
  use_com = True, # apply on centers of mass instead of CA
  colabfold = 1, # PAE format (EBI AF=0, Colabfold=1&2)
  k_harmonic = 700., # Restraint force constant
)

components.add(name=args.fiber)
components.add(name=args.peptide, nmol=105, restraint = False)


components.write(path,name='components.yaml')

