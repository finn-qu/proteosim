from proteosim.mass_spectra_simulation import (
    calculate_mol_mass,
    calculate_mol_mass_collection,
    calculate_mz_collection,
    fragment_peptide,
    amino_acid_mass_dalton
)

def test_calculate_mol_mass():
    aa_mass_dict = {
        "A": 1,
        "B": 2,
        "C": 3,
        "D": 4
        }
    test_peptide1 = "ABCD"
    test_peptide2 = "ABC"
    test_peptide3 = "BD"
    test_peptide4 = "DACB"
    test_peptide5 = "CAB"
    expected_mass1 = {test_peptide1: 10}
    expected_mass2 = {test_peptide2: 6}
    expected_mass3 = {test_peptide3: 6}
    expected_mass4 = {test_peptide4: 10}
    expected_mass5 = {test_peptide5: 6}

    calculated_mass1 = calculate_mol_mass(test_peptide1, aa_mass_dict)
    calculated_mass2 = calculate_mol_mass(test_peptide2, aa_mass_dict)
    calculated_mass3 = calculate_mol_mass(test_peptide3, aa_mass_dict)
    calculated_mass4 = calculate_mol_mass(test_peptide4, aa_mass_dict)
    calculated_mass5 = calculate_mol_mass(test_peptide5, aa_mass_dict)

    assert calculated_mass1 == expected_mass1
    assert calculated_mass2 == expected_mass2
    assert calculated_mass3 == expected_mass3
    assert calculated_mass4 == expected_mass4
    assert calculated_mass5 == expected_mass5


def test_calculate_mol_mass_collection():
    aa_mass_dict = {
            "A": 1,
            "B": 2,
            "C": 3,
            "D": 4
            }
    peptides = ["AA", "AB", "ABC", "BCAD", "BD"]
    expected = {
        "AA": 2,
        "AB": 3,
        "ABC": 6,
        "BCAD": 10,
        "BD": 6
    }
    actual = calculate_mol_mass_collection(peptides, aa_mass_dict)

    assert actual == expected


def test_calculate_mz_collection():
    peptide_mass_map = {
        "AA": 2,
        "AB": 3,
        "BC": 5
        }
    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {
        "AA": (2 + 2 * 1.007) / 2,
        "AB": (3 + 2 * 1.007) / 2,
        "BC": (5 + 2 * 1.007) / 2
    }

    assert actual == expected


def test_fragment_peptide():
    peptide = "PEPT"
    expected = ["P", "PE", "PEP", "T", "PT", "EPT"]
    actual = fragment_peptide(peptide)

    assert set(actual) == set(expected)
    assert "EP" not in actual
    assert "PEPT" not in actual