import numpy as np
import matplotlib.pyplot as plt

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    Sum up the masses of the amino acids in a peptide sequence.

    Parameters
    ----------
    peptide_seq : str
        A string containing the sequence of the peptide.
    amino_acid_mass_dict : dict
        With the amino acids as keys and their molecular masses as values.
    Returns
    -------
    dict
        With the peptide sequence as key and the peptide's molecular mass as value.
    """

    mass = sum(list(map(lambda a: amino_acid_mass_dict[a], peptide_seq)))
    mol_mass = {peptide_seq: mass}
    return mol_mass


def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    For a list of peptides, sum up the masses of their amino acids and return a dictionary of the peptides and their molecular masses.
    
    Parameters
    ----------
    peptides : list
        A list of peptide sequences.
    amino_acid_mass_dict : dict
        With the amino acids as keys and their molecular masses as values.
    Returns
    -------
    dict
        With peptide sequences as keys and molecular masses as values.
    """

    mol_mass_collection = {}
    for peptide in peptides:
        mol_mass_collection.update(calculate_mol_mass(peptide, amino_acid_mass_dict))

    return mol_mass_collection


def calculate_mz_collection(peptide_mass_map, charge, proton_mass=1.007):
    """
    Map the peptide sequences to the m/z ratios obtained after ionization.

    Parameters
    ----------
    peptide_mass_map : dict
        With peptide sequences as keys and molecular masses as values.
    charge : int
        The charge that is added by the ionization?
    proton_mass : float, optional
        The mass of a proton (the default is 1.007).
    Returns
    -------
    dict
        With peptide sequences as keys and m/z ratios as values.
    """

    mz_collection = {}
    for peptide, neutral_mass in peptide_mass_map.items():
        mz_collection[peptide] = (neutral_mass + charge * proton_mass) / charge

    return mz_collection


def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Generate a random intensity for each mz value and create a bar plot of intensities and mz values.

    Parameters
    ----------
    mz_values : list
        A list of mz values.
    random_count_range : tuple
        A tuple containing the upper and lower bound for intensities.
    seed : int
        The seed for random number generation.
    Returns
    -------
    None
    """

    np.random.seed(seed)

    random_lower, random_upper = random_count_range
    spectrum = []
    for mz in mz_values:
        spectrum.append(np.random.randint(random_lower, random_upper))

    plt.bar(mz_values, spectrum, width=1)
    plt.xlabel("m/z value")
    plt.ylabel("intensity")
    plt.show()


def fragment_peptide(peptide):
    """
    Fragment the peptide sequence into the b- and y-ions.

    Parameters
    ----------
    peptide : str
        The sequence of the peptide.
    Returns
    -------
    list
        List containing the sequences of all b- and y-ions
    """
    b_ions = []
    y_ions = []

    for i in range(len(peptide) - 1):
        b_ions.append(peptide[:i+1])
        y_ions.append(peptide[i+1:])
    y_ions.reverse()

    return b_ions + y_ions


amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}