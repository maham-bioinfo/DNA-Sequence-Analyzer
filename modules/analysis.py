from collections import Counter
def sequence_length(sequence):
    """
    Returns the total number of nucleotides.
    """
    sequence = sequence.upper().replace(" ", "").replace("\n", "")
    return len(sequence)
def nucleotide_count(sequence):
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }
def gc_content(sequence):
    
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    length = len(sequence)

    if length == 0:
        return 0

    gc = sequence.count("G") + sequence.count("C")

    gc_percent = (gc / length) * 100

    return round(gc_percent, 2)
def at_content(sequence):
    
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    length = len(sequence)

    if length == 0:
        return 0

    at = sequence.count("A") + sequence.count("T")

    at_percent = (at / length) * 100

    return round(at_percent, 2)
def reverse_complement(sequence):
    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    reverse_comp = ""

    for base in sequence:
        reverse_comp += complement[base]

    return reverse_comp[::-1]
def transcribe(sequence):
    sequence = sequence.upper()
    rna = sequence.replace("T", "U")
    return rna
def translate(sequence):
    rna = sequence.upper().replace("T", "U")

    codon_table = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",

    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",

    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",

    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G"
}

    protein = ""

    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]

        if len(codon) != 3:
            break

        protein += codon_table.get(codon, "X")

    return protein
def amino_acid_frequency(sequence):

    protein = translate(sequence)

    protein = protein.replace("*", "")

    frequency = Counter(protein)

    return dict(frequency)
def find_orf(sequence):
    sequence = sequence.upper()

    start = sequence.find("ATG")

    if start == -1:
        return "No Start Codon Found"

    stop_codons = ["TAA", "TAG", "TGA"]

    for i in range(start + 3, len(sequence), 3):
        codon = sequence[i:i+3]

        if codon in stop_codons:
            return sequence[start:i+3]

    return "No Stop Codon Found"
def read_fasta(file_path):
    sequence = ""

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith(">"):
                continue
            sequence += line.strip()

    return sequence.upper()
def read_multi_fasta(file_path):
    sequences = {}
    header = ""
    sequence = ""

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                if header:
                    sequences[header] = sequence

                header = line[1:]
                sequence = ""
            else:
                sequence += line.upper()

        if header:
            sequences[header] = sequence

    return sequences
def motif_search(sequence, motif):
    sequence = sequence.upper()
    motif = motif.upper()

    positions = []

    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i+len(motif)] == motif:
            positions.append(i + 1)

    return positions
def restriction_sites(sequence):
    sequence = sequence.upper()

    enzymes = {
        "EcoRI": "GAATTC",
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT",
        "NotI": "GCGGCCGC"
    }

    results = {}

    for enzyme, site in enzymes.items():

        positions = []

        for i in range(len(sequence) - len(site) + 1):

            if sequence[i:i+len(site)] == site:
                positions.append(i + 1)

        results[enzyme] = positions

    return results
def save_report(sequence, filename="report.txt"):

    counts = nucleotide_count(sequence)

    with open(filename, "w") as file:

        file.write("DNA Sequence Analyzer Report\n")
        file.write("=" * 40 + "\n\n")

        file.write(f"Sequence:\n{sequence}\n\n")

        file.write(f"Length: {sequence_length(sequence)}\n\n")

        file.write("Nucleotide Counts\n")

        for nucleotide, count in counts.items():
            file.write(f"{nucleotide} : {count}\n")

        file.write(f"\nGC Content : {gc_content(sequence)}%\n")
        file.write(f"AT Content : {at_content(sequence)}%\n\n")

        file.write(f"Reverse Complement:\n{reverse_complement(sequence)}\n\n")

        file.write(f"RNA:\n{transcribe(sequence)}\n\n")

        file.write(f"Protein:\n{translate(sequence)}\n\n")

        file.write(f"ORF:\n{find_orf(sequence)}\n")
def colored_sequence(sequence):

    colors = {
        "A": "green",
        "T": "red",
        "G": "orange",
        "C": "blue"
    }

    html = ""

    for base in sequence.upper():
        color = colors.get(base, "black")
        html += f"<span style='color:{color}; font-weight:bold; font-size:20px'>{base}</span>"

    return html
def codon_frequency(sequence):

    sequence = sequence.upper()

    codons = {}

    for i in range(0, len(sequence) - 2, 3):

        codon = sequence[i:i+3]

        if len(codon) == 3:

            codons[codon] = codons.get(codon, 0) + 1

    return codons