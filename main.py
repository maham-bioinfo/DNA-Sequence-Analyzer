from modules.validation import validate_sequence
from modules.analysis import (
    sequence_length,
    nucleotide_count,
    gc_content,
    at_content,
    reverse_complement,
    transcribe,
    translate,
    find_orf,
    read_fasta,
    read_multi_fasta,
    motif_search,
    restriction_sites,
    save_report
)


def display_results(sequence):
    if not validate_sequence(sequence):
        print("Invalid DNA Sequence")
        return

    print("\nValid DNA Sequence")
    print("Sequence Length:", sequence_length(sequence))

    counts = nucleotide_count(sequence)

    print("\nNucleotide Counts")
    for nucleotide, count in counts.items():
        print(f"{nucleotide} : {count}")

    print(f"\nGC Content : {gc_content(sequence)}%")
    print(f"AT Content : {at_content(sequence)}%")

    print("\nReverse Complement:")
    print(reverse_complement(sequence))

    print("\nRNA Transcription:")
    print(transcribe(sequence))

    print("\nProtein Translation:")
    print(translate(sequence))

    print("\nOpen Reading Frame (ORF):")
    print(find_orf(sequence))

    motif = input("\nEnter motif to search: ")

    positions = motif_search(sequence, motif)

    if positions:
        print("Motif Found at Positions:", positions)
    else:
        print("Motif Not Found")
    print("\nRestriction Enzyme Analysis")

sites = restriction_sites(sequence)

for enzyme, positions in sites.items():

    if positions:
        print(f"{enzyme}: Found at {positions}")
    else:
        print(f"{enzyme}: Not Found")
choice = input("\nDo you want to save the report? (y/n): ")

if choice.lower() == "y":

    filename = input("Enter file name: ")

    save_report(sequence, filename)

    print("Report Saved Successfully!")


print("=" * 40)
print("DNA Sequence Analyzer")
print("=" * 40)

print("1. Enter DNA Sequence")
print("2. Read FASTA File")
print("3. Read Multi-FASTA File")

choice = input("\nChoose option (1/2/3): ")

if choice == "1":

    sequence = input("Enter DNA Sequence: ")
    display_results(sequence)

elif choice == "2":

    file_path = input("Enter FASTA File Path: ")

    try:
        sequence = read_fasta(file_path)
        display_results(sequence)

    except FileNotFoundError:
        print("File Not Found!")

elif choice == "3":

    file_path = input("Enter Multi-FASTA File Path: ")

    try:
        sequences = read_multi_fasta(file_path)

        for header, sequence in sequences.items():

            print("\n" + "=" * 50)
            print("Sequence:", header)
            print("=" * 50)

            display_results(sequence)

    except FileNotFoundError:
        print("File Not Found!")

else:
    print("Invalid Choice")
    