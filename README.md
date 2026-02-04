# aSyn-fibrils-MD
Molecular dynamics simulations of $\alpha$-synuclein fibrils and natural peptides using the [CALVADOS model](https://github.com/KULL-Centre/CALVADOS) (Tesei et al., 2021). This repository contains specific data and coding scripts to run CALVADOS with amyloid fibrils and peptides. 

## Documentation
Topology and trajectory files of modelled systems to be found in CORA.RDR repository. CORA repository contains both topology (.pdb) and trajectory (.dcd) files obtained for all simulations.

### File system directories

- Folder name: `./aggrescan`
  - Description: input and output files for the prediction of aggregation-prone regions with AGGRESCAN and Aggrescan4D.

- Folder name: `./all_peptides`
  - Description: csv files with predicted contacts and properties for aSynPEP-DB peptides and negative entries.

- Folder name: `./analysis`
  - Description: Python scripts to analyze trajectories and topology files.

- Folder name: `./input`
  - Description: original structures and sequences of peptides and fibrils.

- Folder name: `./contact_arrays`
  - Description: data arrays for total contacts of simulated systems (.npy).

- Folder name: `./controls`
  - Description: peptide contacts for control simulations.

- Folder name: `./running_CALVADOS`
  - Description: input files to prepare and generate specific simulations of peptide-fibrils. Please refer to original [CALVADOS repository](https://github.com/KULL-Centre/CALVADOS) for general guidelines on running the model.



## References

### CALVADOS model

For general usage of CALVADOS please refer to the original works/repositories of the authors:

- G. Tesei, T. K. Schulze, R. Crehuet, K. Lindorff-Larsen. Accurate model of liquid-liquid phase behavior of intrinsically disordered proteins from optimization of single-chain properties. PNAS (2021), 118(44):e2111696118. DOI: [10.1073/pnas.2111696118](https://www.pnas.org/doi/full/10.1073/pnas.2111696118).

- G. Tesei, K. Lindorff-Larsen. Improved predictions of phase behaviour of intrinsically disordered proteins by tuning the interaction range. Open Research Europe (2022), 2(94). DOI: [10.12688/openreseurope.14967.2](https://open-research-europe.ec.europa.eu/articles/2-94/v2).

- F. Cao, S. von Bülow, G. Tesei, K. Lindorff-Larsen. A coarse-grained model for disordered and multi-domain proteins. Protein Science (2024), 33(11):e5172. DOI: [10.1002/pro.5172](https://onlinelibrary.wiley.com/doi/10.1002/pro.5172).

- S. von Bülow*, Y. Yasuda#, F. Cao#, T. K. Schulze#, A. I. Trolle#, A. S. Rauh#, R. Crehuet#, K. Lindorff-Larsen*, G. Tesei* (# equal contribution) Software package for simulations using the coarse-grained CALVADOS model, arXiv 2025. [https://doi.org/10.48550/arXiv.2504.10408](https://doi.org/10.48550/arXiv.2504.10408)


### Peptides and $\alpha$-synuclein

For referencing research related to peptides and $\alpha$-synuclein, please cite: 

- Pintado-Grima C, Bárcenas O, Iglesias V, Santos J, Manglano-Artuñedo Z, Pallarès I, Burdukiewicz M, Ventura S. aSynPEP-DB: a database of biogenic peptides for inhibiting α-synuclein aggregation. Database, Volume 2023, 2023. DOI: [10.1093/database/baad084](https://academic.oup.com/database/article/doi/10.1093/database/baad084/7451591?login=false).

- Santos J, Gracia P, Navarro S, Peña-Díaz S, Pujols J, Cremades N, Pallarès I, Ventura S. α-Helical peptidic scaffolds to target α-synuclein toxic species with nanomolar affinity. Nat Commun. 2021 Jun 18;12(1):3752. DOI: [10.1038/s41467-021-24039-2](https://www.nature.com/articles/s41467-021-24039-2).

- Pintado-Grima C, Ventura S. The role of amphipathic and cationic helical peptides in Parkinson's disease. Protein Sci. 2025. DOI: [10.1002/pro.70020](https://onlinelibrary.wiley.com/doi/10.1002/pro.70020).

- Santos J, Pallarès I, Ventura S. Is a cure for Parkinson’s disease hiding inside us? Trends Biochem Sci. 2022 Aug;47(8):641-644. DOI: [10.1016/j.tibs.2022.02.001](https://www.sciencedirect.com/science/article/pii/S0968000422000251).

