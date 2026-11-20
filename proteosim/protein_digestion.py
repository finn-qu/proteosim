import re

enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
}

def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Split the protein sequence at the cleavage pattern and return all the obtained peptides that have the correct length.

    Parameters
    ----------
    protein_seq : str
            A string containing the amino acid sequence.
    cleave_pattern: str
            Regular expression as the cleavage pattern.
    min_pep_len : int, optional
        Minimum number of amino acids a peptide needs to have to be returned (the default is 5).
    max_pep_len : int, optional
        Maximum number of amino acids a peptide needs to have to be returned (the default is 30).
    Returns
    -------
    list of str
        List containing the substrings corresponding to the peptides obtained after digestion.
    """

    peptides = re.split(cleave_pattern, protein_seq)
    filtered_peptides = list(filter(lambda a : len(a) >= min_pep_len and len(a) <= max_pep_len, peptides))
    return filtered_peptides

def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    
    """
    Split all the protein sequences at the cleavage pattern and return a dictionary mapping the protein ids to all the 
    corresponding peptides obtained.
    
    Parameters
    ----------
    protein_map : dict
            A dictionary mapping protein ids to sequences.
    cleave_pattern: str
        Regular expression as the cleavage pattern.
    min_pep_len : int, optional
        Minimum number of amino acids a peptide needs to have to be returned (the default is 5).
    max_pep_len : int, optional
        Maximum number of amino acids a peptide needs to have to be returned (the default is 30).
    Returns
    -------
    dict
        Dictionary mapping protein ids to peptide sequences obtained after digestion.
    """
    
    peptide_collection = {}
    for id, seq in protein_map.items():
        peptide_collection[id] = digest_protein_sequence(seq, cleave_pattern, min_pep_len, max_pep_len)

    return peptide_collection

def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute the percentage of amino acids represented by the peptides after digestion. Compare length of whole protein to combined lengths
    of all peptides. Works only if the peptides are correct and cut only by one protease.

    Parameters
    ----------
    protein_seq : str
        The protein's amino acid sequence.
    peptides : list
        A list containing the peptides sequences as strings.
    Returns
    -------
    float
        Percentage of amino acids represented by the peptides after digestion.
    """

    # the coverage depends on how many peptides are filtered out due to being too short or too long
    # so if the coverage is low, it is because fragments are "missing" 
    # the combined length of all peptides tells us how many amino acids are represented by the peptides after digestion

    # calculate the length of each fragment
    peptide_lengths = map(lambda a: len(a), peptides)

    # sum up the length of all fragments and divide it by the length of the whole sequence to get the coverage
    coverage = sum(peptide_lengths) / len(protein_seq) * 100

    return coverage