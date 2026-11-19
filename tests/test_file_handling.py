from pathlib import Path
from proteosim.file_handling import read_fasta

def test_fasta_reader():
    tmp_fasta_path = Path(__file__).parent / '../data/dummy_proteins.fasta'
    protein_map = read_fasta(tmp_fasta_path)
    
    # Replace the strings with your fasta content
    # which you expect to be now available as a dictionary
    assert protein_map["ID1"] == "BLABLABLA"
    assert protein_map["ID2"] == "TRALALALALALA"

test_fasta_reader()