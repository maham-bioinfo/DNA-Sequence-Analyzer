VALID_BASES = {"A", "T", "G", "C"}

def validate_sequence(sequence):
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    for base in sequence:
        if base not in VALID_BASES:
            return False

    return True