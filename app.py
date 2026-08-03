import streamlit as st
import hashlib
import pandas as pd
from modules.validation import validate_sequence
from modules.analysis import read_fasta
from ui import analyze_sequence

st.set_page_config(
    page_title="DNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 DNA Sequence Analyzer")

st.sidebar.title("📂 Input Method")

option = st.sidebar.radio(
    "Choose Input",
    (
        "Paste DNA Sequence",
        "Upload Single FASTA",
        "Upload Multiple FASTA"
    )
)

sequence = ""

uploaded_files = None

# ---------------------------------
# Paste DNA
# ---------------------------------

if option == "Paste DNA Sequence":

    sequence = st.text_area(
        "Paste DNA Sequence",
        height=180,
        placeholder="ATGCGATCGATCG..."
    )

# ---------------------------------
# Upload Single FASTA
# ---------------------------------

elif option == "Upload Single FASTA":

    uploaded_file = st.file_uploader(
        "Upload FASTA File",
        type=["fasta", "fa", "txt"]
    )

    if uploaded_file:

        with open("temp.fasta", "wb") as f:
            f.write(uploaded_file.getbuffer())

        sequence = read_fasta("temp.fasta")

        st.success("✅ FASTA Loaded Successfully")

# ---------------------------------
# Upload Multiple FASTA
# ---------------------------------

else:

    uploaded_files = st.file_uploader(
        "Upload Multiple FASTA Files",
        type=["fasta", "fa", "txt"],
        accept_multiple_files=True
    )

motif = st.text_input(
    "Motif (Optional)",
    placeholder="Example: ATG"
)

analyze = st.button("Analyze")
# ============================================================
# MULTIPLE FASTA ANALYSIS
# ============================================================

if analyze and option == "Upload Multiple FASTA":

    if not uploaded_files:

        st.warning("Please upload one or more FASTA files.")
        st.stop()

    st.header("🧬 Multi FASTA Analysis")

    summary = []

    for uploaded_file in uploaded_files:

        uploaded_file.seek(0)

        text = uploaded_file.read().decode("utf-8")

        dna = ""

        for line in text.splitlines():

            if not line.startswith(">"):
                dna += line.strip()

        dna = dna.upper()

        counts = {
            "A": dna.count("A"),
            "T": dna.count("T"),
            "G": dna.count("G"),
            "C": dna.count("C")
        }

        summary.append({
            "File": uploaded_file.name,
            "Length": len(dna),
            "GC %": round((counts["G"] + counts["C"]) / len(dna) * 100, 2),
            "AT %": round((counts["A"] + counts["T"]) / len(dna) * 100, 2),
            "A": counts["A"],
            "T": counts["T"],
            "G": counts["G"],
            "C": counts["C"]
        })

        st.write(f"Sequence length: {len(dna)} bp")

        if not validate_sequence(dna):

            st.error(f"❌ {uploaded_file.name} contains invalid DNA sequence.")
            continue

        st.divider()

        st.subheader(f"📄 {uploaded_file.name}")

        unique_key = hashlib.md5(dna.encode()).hexdigest()

        analyze_sequence(
            dna=dna,
            motif=motif,
            unique_key=unique_key
        )
    st.header("📊 FASTA Statistics Summary")

    summary_df = pd.DataFrame(summary)

    st.dataframe(summary_df, width="stretch")
         
    st.download_button(
        "📥 Download Summary CSV",
        summary_df.to_csv(index=False),
        file_name="FASTA_Summary.csv",
        mime="text/csv"
    )
    st.stop()
# ============================================================
# SINGLE FASTA / PASTED DNA ANALYSIS
# ============================================================

if analyze and option != "Upload Multiple FASTA":

    sequence = sequence.upper().strip()

    if sequence == "":

        st.warning("Please enter a DNA sequence or upload a FASTA file.")
        st.stop()

    if not validate_sequence(sequence):

        st.error("❌ Invalid DNA Sequence")
        st.stop()

    st.success("✅ Valid DNA Sequence")

    analyze_sequence(
        dna=sequence,
        motif=motif,
        unique_key="single_sequence"
    )