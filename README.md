# Proteosim Course Package
This repository allows you to simulate a proteomics pipeline. It includes functions for protein digestion by proteases, Liquid Chromatography and Mass Spectrometry. The code is under active development.

## Functionality
### Protein digestion
- simulate a proteomics pipeline
- load & parse a FASTA file
- digest the proteins using different proteases & filters the peptide lengths
- calculate the coverage of a protein by the peptides obtained after digestion by different enzymes
- compare the coverage of all protein-enzyme-combinations & creates a boxplot for visualization


### Liquid chromatography
- predict the LC retention times of the peptides based on `achrom.RCs_guo_ph7_0`
- create a histogram of the retention times
- filter the peptides based on retention times


### Mass Spectra Simulation
- calculate the molecular mass of a peptide or a set of peptides based on the molecular masses of the amino acids
- calculate the m/z values of the peptides after ionization
- plot a spectrum of m/z values based on randomly generated intensities
- fragment a peptide into b- and y-ions

## Requirements & Installation
### Requirements
- numpy
- pandas
- matplotlib
- pyteomics


### Installation
- idk


## Usage