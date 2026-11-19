def read_fasta(filepath):
    """
    Load a FASTA file and create a dictionary that maps each protein id to the amino acid sequence

    Parameters
    filepath : str
        The path of the FASTA file.
    Returns
    -------
    dict
        With protein ids as keys and sequences as values.
    """
    
    # creates an empty dictionary that will store the identifiers and sequences, create an empty list that will store the sequence
    protein_map = {}
    current_id = None
    current_sequence = []

    # open the file
    with open(filepath, 'r') as fasta_handle:

        # go through the file line by line
        for line in fasta_handle:

            # remove spaces
            stripped = line.strip()

            # if a line is empty, skip the next steps
            if not stripped:
                continue

            # if a line is a header and has an id, add the id to the dictionary and join the previous sequence together to create a string, then add it to the dictionary
            if stripped.startswith('>'):
                if current_id is not None:
                    protein_map[current_id] = ''.join(current_sequence)
                    current_sequence = []

                # get the id from the current header
                current_id = stripped.split('|')[1]

            # if a line is not a header, add it to the list for the dictionary
            else:
                current_sequence.append(stripped)

    # join the entries in the sequence list to create a string
    if current_id is not None:
        protein_map[current_id] = ''.join(current_sequence)

    return protein_map