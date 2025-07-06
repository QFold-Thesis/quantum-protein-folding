def aa_to_hp(sequence: str) -> str:
    """
    Converts an amino acid sequence to a hydrophobic/polar (HP) sequence.
    """
    # List of hydrophobic amino acids according to the most commonly used HP convention
    hydrophobic = {'A', 'V', 'L', 'I', 'M', 'F', 'W', 'Y', 'C'}
    # The rest are considered polar
    hp_sequence = ""
    for aa in sequence.upper():
        if aa in hydrophobic:
            hp_sequence += "H"
        else:
            hp_sequence += "P"
    return hp_sequence