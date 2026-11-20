# Proteosim Course Package
This repository allows you to simulate a proteomics pipeline. It includes functions for protein digestion by proteases, Liquid Chromatography and Mass Spectrometry. The code is under active development.

## Usage
- input: FASTA file (proteins + sequences)
- processing steps:
  - simulation of digestion
  - Liquid Chromatography (LC) simulation
  - Mass Spectrometry (MS) simulation
- outputs: 
  - m/z values of peptides & peptide fragments
  - histogram of LC retention times
  - MS1 & MS2 spectra

## Requirements & Installation
### Requirements
- numpy
- pandas
- matplotlib
- pyteomics


### Installation
`pip install -r requirements.txt`  
`pip install .`


## Functions
### Protein digestion
- read_fasta: load & parse a FASTA file
- digest_protein_sequence: digest a protein using different proteases & filter the peptide length
- digest_protein_collection: digest a dictionary of proteins using different proteases & filter the peptide lengths
- compute_sequence_coverage: calculate the coverage of a protein by the peptides obtained after digestion
- MISSING: compare the coverage of all protein-enzyme-combinations & creates a boxplot for visualization


### Liquid chromatography
- predict_lc_retention_times: predict the LC retention times of a list of peptides based on `achrom.RCs_guo_ph7_0`
- plot_retention_time: create a histogram of the retention times
- select_retention_time_window: filter the peptides based on retention times


### Mass Spectra Simulation
- calculate_mol_mass: calculate the molecular mass of a peptide based on the molecular masses of the amino acids
- calculate_mol_mass_collection: calculate the molecular mass of a set of peptides based on the molecular masses of the amino acids
- calculate_mz_collection: calculate the m/z values of the peptides after ionization
- plot_spectrum: plot a spectrum of m/z values based on randomly generated intensities
- fragment_peptide: fragment a peptide into b- and y-ions


## Example of Usage
- for the quickest way to try the package, see tutorials/ms_experiment_final.ipynb for a step-by-step usage of the package