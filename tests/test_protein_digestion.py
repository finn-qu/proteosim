import pandas as pd
from proteosim.protein_digestion import (
    digest_protein_sequence, 
    digest_protein_collection, 
    compute_sequence_coverage, 
    compare_sequence_coverages,
    enzyme_cleavage_patterns
    )

def test_digest_protein_sequence():
    dummy_seq_1 = "ABKABRABKPA"

    assert digest_protein_sequence(dummy_seq_1, enzyme_cleavage_patterns["LysC"], min_pep_len=0, max_pep_len=30) == ["ABK","ABRABK","PA"]
    assert digest_protein_sequence(dummy_seq_1, enzyme_cleavage_patterns["LysC"], min_pep_len=3, max_pep_len=30) == ["ABK", "ABRABK"]
    assert digest_protein_sequence(dummy_seq_1, enzyme_cleavage_patterns["LysC"], min_pep_len=4, max_pep_len=30) == ["ABRABK"]
    assert digest_protein_sequence(dummy_seq_1, enzyme_cleavage_patterns["Trypsin"], min_pep_len=0, max_pep_len=30) == ["ABK", "ABR", "ABKPA"]
    assert digest_protein_sequence(dummy_seq_1, enzyme_cleavage_patterns["Trypsin"], min_pep_len=4, max_pep_len=30) == ["ABKPA"]


def test_digest_protein_collection():
    dummy_proteins = {
        "ID3": "ABKABRABKPA",
        "ID4": "TRAKKARPAKPA"
    }

    assert digest_protein_collection(dummy_proteins, enzyme_cleavage_patterns["Trypsin"], min_pep_len=0, max_pep_len=30) == {
        "ID3" : ["ABK","ABR","ABKPA"],
        "ID4" : ["TR","AK", "K","ARPAKPA"]
        }
    assert digest_protein_collection(dummy_proteins, enzyme_cleavage_patterns["Trypsin"], min_pep_len=2,max_pep_len=30) == {
        "ID3" : ["ABK","ABR","ABKPA"],
        "ID4" : ["TR","AK","ARPAKPA"]
        }
    assert digest_protein_collection(dummy_proteins, enzyme_cleavage_patterns["Trypsin"], min_pep_len=3,max_pep_len=30) == {
        "ID3" : ["ABK","ABR","ABKPA"],
        "ID4" : ["ARPAKPA"]
        }


def test_compute_sequence_coverage():
    dummy_prot_seq = "ABKABRABKP"
    dummy_peps_1 = ["ABK", "ABR", "ABKP"]
    dummy_peps_2 = ["AB"]
    dummy_peps_3 = ["ABK", "RAB", "KP"]

    coverage_1 = compute_sequence_coverage(dummy_prot_seq, dummy_peps_1)
    coverage_2 = compute_sequence_coverage(dummy_prot_seq, dummy_peps_2)
    coverage_3 = compute_sequence_coverage(dummy_prot_seq, dummy_peps_3)
    assert coverage_1 == 100
    assert coverage_2 == 20
    assert coverage_3 == 80

def test_compare_sequence_coverages():
    dummy_protein_map = {
        "ID1": "AKWAWAWAWA",
        "ID2": "ARPAKWAWAW"
    }
    dummy_peptides_trypsin = {
        "ID1": ["AK", "WAWAWAWA"],
        "ID2": ["ARPAK", "WAWAW"]
    }
    expected = pd.DataFrame({
        "ID1": 80,
        "ID2": 100
    })

    #assert compare_sequence_coverages(dummy_protein_map, enzyme_cleavage_patterns) == expected
    #print(expected.eq(compare_sequence_coverages(dummy_protein_map, enzyme_cleavage_patterns)))
    print(expected)
    print(compare_sequence_coverages(dummy_protein_map, enzyme_cleavage_patterns))

