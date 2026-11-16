import os
import subprocess
import tempfile
import shutil
from Bio import AlignIO
from Bio.Align import PairwiseAligner
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
def gc_sliding_window(sequence, window_size=100):

    sequence = sequence.upper()

    positions = []
    gc_values = []

    for i in range(0, len(sequence) - window_size + 1):

        window = sequence[i:i + window_size]

        gc = (
            window.count("G") +
            window.count("C")
        ) / window_size * 100

        positions.append(i + 1)
        gc_values.append(round(gc, 2))

    return positions, gc_values
def dna_quality_score(sequence):

    sequence = sequence.upper()

    score = 0
    report = []

    # Valid DNA
    valid = set(sequence).issubset({"A", "T", "G", "C"})

    if valid:
        score += 20
        report.append(("✅ Valid DNA", True))
    else:
        report.append(("❌ Invalid DNA", False))

    # GC %

    gc = gc_content(sequence)

    if 40 <= gc <= 60:
        score += 20
        report.append(("✅ Normal GC Content", True))
    else:
        report.append(("⚠ Unusual GC Content", False))

    # Start codon

    if "ATG" in sequence:
        score += 20
        report.append(("✅ Start Codon Found", True))
    else:
        report.append(("❌ Start Codon Missing", False))

    # Stop codon

    if any(x in sequence for x in ["TAA","TAG","TGA"]):
        score += 20
        report.append(("✅ Stop Codon Found", True))
    else:
        report.append(("❌ Stop Codon Missing", False))

    # ORF

    if "No" not in find_orf(sequence):
        score += 20
        report.append(("✅ ORF Found", True))
    else:
        report.append(("❌ ORF Missing", False))

    return score, report
def genome_features(sequence, motif=""):

    sequence = sequence.upper()

    features = []

    # ORF
    start = sequence.find("ATG")

    if start != -1:

        stop_codons = ["TAA", "TAG", "TGA"]

        for i in range(start + 3, len(sequence), 3):

            codon = sequence[i:i+3]

            if codon in stop_codons:

                features.append({
                    "Feature": "ORF",
                    "Start": start + 1,
                    "End": i + 3
                })

                break

    # Restriction Sites

    enzymes = restriction_sites(sequence)

    for enzyme, positions in enzymes.items():

        for pos in positions:

            features.append({
                "Feature": enzyme,
                "Start": pos,
                "End": pos + 5
            })

    # Motif

    if motif:

        positions = motif_search(sequence, motif)

        for pos in positions:

            features.append({
                "Feature": "Motif",
                "Start": pos,
                "End": pos + len(motif) - 1
            })

    return features
def genome_map_data(sequence, motif=""):

    sequence = sequence.upper()

    data = []

    # ORF
    start = sequence.find("ATG")

    if start != -1:

        stop_codons = ["TAA", "TAG", "TGA"]

        for i in range(start + 3, len(sequence), 3):

            codon = sequence[i:i+3]

            if codon in stop_codons:

                data.append({
                    "type": "ORF",
                    "start": start + 1,
                    "end": i + 3
                })

                break

    # Restriction Sites

    enzymes = restriction_sites(sequence)

    for enzyme, positions in enzymes.items():

        for p in positions:

            data.append({
                "type": enzyme,
                "start": p,
                "end": p
            })

    # Motif

    if motif:

        for p in motif_search(sequence, motif):

            data.append({
                "type": "Motif",
                "start": p,
                "end": p + len(motif)
            })

    return data
def genome_dataframe(sequence, motif=""):

    import pandas as pd

    features = genome_map_data(sequence, motif)

    rows = []

    tracks = {
        "ORF": 1,
        "Motif": 2,
        "EcoRI": 3,
        "BamHI": 4,
        "HindIII": 5,
        "NotI": 6
    }

    colors = {
        "ORF": "green",
        "Motif": "orange",
        "EcoRI": "red",
        "BamHI": "blue",
        "HindIII": "purple",
        "NotI": "brown"
    }

    for f in features:

        rows.append({

            "Feature": f["type"],

            "Start": f["start"],

            "End": f["end"],

            "Track": tracks.get(f["type"],7),

            "Color": colors.get(f["type"],"black")

        })

    return pd.DataFrame(rows)
def pairwise_alignment(
    seq1,
    seq2,
    mode="global",
    match=2,
    mismatch=-1,
    gap_open=-2,
    gap_extend=-0.5
):

    aligner = PairwiseAligner()

    aligner.mode = mode
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend

    alignments = aligner.align(seq1.upper(), seq2.upper())

    return alignments[0]
def run_mafft_msa(sequences):
    """
    Run MAFFT Multiple Sequence Alignment.

    Parameters
    ----------
    sequences : dict
        Dictionary containing sequence names and DNA sequences.

    Returns
    -------
    MultipleSeqAlignment
        MAFFT aligned sequences.
    """

    if len(sequences) < 2:
        raise ValueError(
            "At least 2 sequences are required for MSA."
        )

    # ---------------------------------
    # Find MAFFT executable
    # ---------------------------------

    mafft_path = shutil.which("mafft")

    if mafft_path is None:

        possible_paths = [
            r"C:\Program Files\MAFFT\mafft.bat",
            r"C:\Program Files\mafft-win\mafft.bat",
            r"C:\Program Files\MAFFT\mafft-win.bat"
        ]

        for path in possible_paths:

            if os.path.exists(path):

                mafft_path = path
                break

    if mafft_path is None:

        raise FileNotFoundError(
            "MAFFT was not found. "
            "Please install MAFFT or add it to PATH."
        )

    # ---------------------------------
    # Create temporary FASTA file
    # ---------------------------------

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".fasta",
        delete=False,
        encoding="utf-8"
    ) as fasta_file:

        fasta_path = fasta_file.name

        for name, sequence in sequences.items():

            fasta_file.write(
                f">{name}\n"
            )

            fasta_file.write(
                f"{sequence}\n"
            )

    try:

        # ---------------------------------
        # Run MAFFT
        # ---------------------------------

        result = subprocess.run(
            [
                mafft_path,
                "--auto",
                fasta_path
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace"
        )

        if result.returncode != 0:

            raise RuntimeError(
                "MAFFT failed:\n"
                + result.stderr
            )

        # ---------------------------------
        # Parse MAFFT output
        # ---------------------------------

        alignment = AlignIO.read(
            result.stdout.splitlines(),
            "fasta"
        )

        return alignment

    finally:

        # ---------------------------------
        # Remove temporary file
        # ---------------------------------

        if os.path.exists(fasta_path):

            os.remove(fasta_path)