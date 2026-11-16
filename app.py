import streamlit as st
import hashlib
import pandas as pd
from modules.validation import validate_sequence
from modules.analysis import (
    read_fasta,
    pairwise_alignment,
    run_mafft_msa
)
from ui import analyze_sequence, show_pairwise_alignment, show_msa_alignment

st.set_page_config(
    page_title="GenScope | DNA Sequence Analysis",
    page_icon="🧬",
    layout="wide"
)
# ============================================================
# NUCLEOSCOPE — PROFESSIONAL BRANDING & THEME
# ============================================================

st.markdown("""
<style>

/* Main application background */
[data-testid="stAppViewContainer"] {
    background: linear-gradient(
        120deg,
        #C4D4F0 0%,
        #CFC1E8 50%,
        #BBDCCF 100%
    );
}

/* Keep Streamlit's main area transparent so the gradient shows */
[data-testid="stMain"],
[data-testid="stHeader"] {
    background: transparent;
}

/* Sidebar gradient */
[data-testid="stSidebar"] {
    background: linear-gradient(
        180deg,
        #AFC2E2 0%,
        #BEB0D9 100%
    );
    border-right: 1px solid rgba(48, 66, 103, 0.12);
}

/* Sidebar brand panel */
.brand-panel {
    text-align: center;
    padding: 18px 10px 24px 10px;
    margin: 4px 0 22px 0;
    background: linear-gradient(
        145deg,
        rgba(255, 255, 255, 0.48),
        rgba(255, 255, 255, 0.20)
    );
    border: 1px solid rgba(255, 255, 255, 0.55);
    border-radius: 18px;
    box-shadow: 0 8px 22px rgba(47, 63, 99, 0.08);
}

.brand-icon {
    font-size: 42px;
    line-height: 1.3;
    margin-bottom: 5px;
}

.brand-name {
    color: #243858;
    font-size: 23px;
    font-weight: 800;
    letter-spacing: 0.2px;
}

.brand-tagline {
    color: #465D80;
    font-size: 11px;
    font-weight: 600;
    letter-spacing: 1.1px;
    text-transform: uppercase;
    margin-top: 7px;
}

.brand-version {
    display: inline-block;
    margin-top: 12px;
    padding: 4px 10px;
    border-radius: 20px;
    background: rgba(36, 56, 88, 0.10);
    color: #30496F;
    font-size: 10px;
    font-weight: 700;
}

/* Sidebar section label */
.sidebar-section-label {
    color: #344B70;
    font-size: 11px;
    font-weight: 800;
    text-transform: uppercase;
    letter-spacing: 1.4px;
    margin: 4px 0 8px 2px;
}

/* Sidebar radio labels */
[data-testid="stSidebar"] [data-testid="stRadio"] label {
    color: #263D61;
    font-weight: 550;
}

/* Main welcome banner */
.dashboard-banner {
    background: linear-gradient(
        120deg,
        #344F79 0%,
        #5A527F 55%,
        #376B69 100%
    );
    border: 1px solid rgba(255, 255, 255, 0.35);
    border-radius: 20px;
    padding: 25px 30px;
    margin: 10px 0 28px 0;
    box-shadow: 0 12px 28px rgba(43, 57, 91, 0.16);
}

.banner-title {
    font-size: 1.55rem;
    font-weight: 800;
    color: #FFFFFF;
    letter-spacing: -0.3px;
}

.banner-text {
    margin-top: 7px;
    color: #E2E9F5;
    font-size: 0.96rem;
    line-height: 1.55;
}

.banner-caption {
    margin-top: 14px;
    color: #D4E2F0;
    font-size: 0.75rem;
    font-weight: 650;
    letter-spacing: 1.2px;
    text-transform: uppercase;
}

.banner-icon {
    float: right;
    font-size: 44px;
    margin-left: 18px;
}

/* Buttons */
.stButton > button,
.stDownloadButton > button {
    border-radius: 10px;
    border: 1px solid #526D98;
    font-weight: 650;
    transition: all 0.2s ease;
}

.stButton > button:hover,
.stDownloadButton > button:hover {
    border-color: #344F79;
    color: #263D61;
    box-shadow: 0 4px 12px rgba(43, 57, 91, 0.12);
}

/* Input fields */
[data-testid="stTextArea"] textarea,
[data-testid="stTextInput"] input {
    border-radius: 12px;
}

/* Subtle divider */
[data-testid="stSidebar"] hr {
    border-color: rgba(48, 66, 103, 0.16);
}

</style>
""", unsafe_allow_html=True)


# Sidebar branding
st.sidebar.markdown(
    """
    <div class="brand-panel">
        <div class="brand-icon">🧬</div>
        <div class="brand-name">GenScope</div>
        <div class="brand-tagline">Explore. Analyze. Understand.</div>
        <div class="brand-version">BIOINFORMATICS PLATFORM · v1.0</div>
    </div>
    """,
    unsafe_allow_html=True
)

st.sidebar.markdown(
    '<div class="sidebar-section-label">Analysis Workspace</div>',
    unsafe_allow_html=True
)

option = st.sidebar.radio(
    "Choose Input",
    (
        "Paste DNA Sequence",
        "Upload Single FASTA",
        "Upload Multiple FASTA",
        "Pairwise Sequence Alignment",
        "Multiple Sequence Alignment"
    )
)
# ---------------------------------
# DASHBOARD BANNER
# ---------------------------------

st.markdown(
    """
    <div class="dashboard-banner">
        <div class="banner-icon">🧬</div>
        <div class="banner-title">Welcome to GenScope</div>
        <div class="banner-text">
            An interactive environment to explore, analyze, compare,
            and visualize nucleotide sequences using bioinformatics tools.
        </div>
        <div class="banner-caption">DNA Sequence Analysis · FASTA Processing · Sequence Alignment</div>
    </div>
    """,
    unsafe_allow_html=True
)

sequence = ""
sequence2 = ""
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

elif option == "Upload Multiple FASTA":

    uploaded_files = st.file_uploader(
        "Upload Multiple FASTA Files",
        type=["fasta", "fa", "txt"],
        accept_multiple_files=True
    )
# ---------------------------------
# Multiple Sequence Alignment
# ---------------------------------

elif option == "Multiple Sequence Alignment":

    st.subheader("🧬 Multiple Sequence Alignment")

    msa_input_method = st.radio(
        "Choose Input Method",
        [
            "Paste FASTA",
            "Upload FASTA File"
        ],
        horizontal=True
    )

    msa_sequences = {}
    # ---------------------------------
    # Paste FASTA
    # ---------------------------------

    if msa_input_method == "Paste FASTA":

        fasta_text = st.text_area(
            "Paste FASTA sequences",
            height=250,
            placeholder=(
                ">Sequence_1\n"
                "ATGCGATCGATCG\n"
                ">Sequence_2\n"
                "ATGCAATCGATCG\n"
                ">Sequence_3\n"
                "ATGCGATGGATCG"
            )
        )

        if fasta_text.strip():

            current_header = None
            current_sequence = []

            for line in fasta_text.splitlines():

                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):

                    if current_header is not None:
                        msa_sequences[current_header] = "".join(
                            current_sequence
                        ).upper()

                    current_header = line[1:].strip()
                    current_sequence = []

                else:

                    current_sequence.append(line)

            if current_header is not None:

                msa_sequences[current_header] = "".join(
                    current_sequence
                ).upper()
    # ---------------------------------
    # Upload FASTA
    # ---------------------------------

    else:

        uploaded_msa = st.file_uploader(
            "Upload FASTA Files",
            type=["fasta", "fa", "fas", "txt"],
            accept_multiple_files=True,
            key="msa_file_uploader"
        )

        if uploaded_msa:

            for uploaded_file in uploaded_msa:

                fasta_text = uploaded_file.read().decode("utf-8")

                current_header = None
                current_sequence = []

                for line in fasta_text.splitlines():

                    line = line.strip()

                    if not line:
                        continue

                    if line.startswith(">"):

                        if current_header is not None:

                            msa_sequences[current_header] = "".join(
                                current_sequence
                            ).upper()

                        current_header = line[1:].strip()
                        current_sequence = []

                    else:

                        current_sequence.append(line)

                if current_header is not None:

                    msa_sequences[current_header] = "".join(
                        current_sequence
                    ).upper()
    # ---------------------------------
    # Validation / Summary
    # ---------------------------------

    if msa_sequences:

        st.success(
            f"✅ {len(msa_sequences)} sequences loaded"
        )

        for name, seq in msa_sequences.items():

            if not validate_sequence(seq):
                st.error(
                    f"❌ Invalid DNA sequence: {name}"
                )
                st.stop()

        if len(msa_sequences) < 2:

            st.warning(
                "⚠️ Multiple Sequence Alignment requires "
                "at least 2 sequences."
            )
            st.stop()

        st.write(
            f"**Sequences:** {len(msa_sequences)}"
        )

        st.write(
            f"**Average length:** "
            f"{sum(len(seq) for seq in msa_sequences.values()) / len(msa_sequences):.1f} bp"
        )

        if st.button(
            "🧬 Run Multiple Sequence Alignment",
            key="run_msa"
        ):

            st.session_state.msa_sequences = msa_sequences
            
            with st.spinner(
                "Running MAFFT Multiple Sequence Alignment..."
            ):

                try:

                    alignment = run_mafft_msa(
                        msa_sequences
                    )

                    st.session_state.msa_alignment = alignment
                    st.session_state.msa_ready = True
                    st.success(
                        "✅ Multiple Sequence Alignment completed successfully!"
                    )

                except FileNotFoundError as e:

                    st.error(
                        f"❌ {e}"
                    )

                except Exception as e:

                    st.error(
                        f"❌ MSA failed: {e}"
                    )
# ---------------------------------
# Pairwise Sequence Alignment
# ---------------------------------

elif option == "Pairwise Sequence Alignment":

    sequence = st.text_area(
        "Sequence 1",
        height=150,
        key="seq1"
    )

    sequence2 = st.text_area(
        "Sequence 2",
        height=150,
        key="seq2"
    )
motif = st.text_input(
    "Motif (Optional)",
    placeholder="Example: ATG"
)

if "analyzed" not in st.session_state:
    st.session_state.analyzed = False

if st.button("Analyze"):
    st.session_state.analyzed = True
    
# Agar input method change ho jaye to reset analysis
if "last_option" not in st.session_state:
    st.session_state.last_option = option

if st.session_state.last_option != option:
    st.session_state.analyzed = False
    st.session_state.last_option = option
# ============================================================
# MSA RESULT
# ============================================================

if (
    option == "Multiple Sequence Alignment"
    and st.session_state.get("msa_ready", False)
):

    show_msa_alignment(
        st.session_state.msa_alignment
    )
# ============================================================
# PAIRWISE ALIGNMENT
# ============================================================

if st.session_state.analyzed and option == "Pairwise Sequence Alignment":

    

    seq1 = sequence.upper().strip()
    seq2 = sequence2.upper().strip()

    if seq1 == "" or seq2 == "":
        st.warning("Please enter both sequences.")
        st.stop()

    if not validate_sequence(seq1):
        st.error("Sequence 1 is invalid.")
        st.stop()

    if not validate_sequence(seq2):
        st.error("Sequence 2 is invalid.")
        st.stop()

    st.session_state.alignment = pairwise_alignment(seq1, seq2)
    

    show_pairwise_alignment(
        st.session_state.alignment
    )
# ============================================================
# MULTIPLE FASTA ANALYSIS
# ============================================================

if st.session_state.analyzed and option == "Upload Multiple FASTA":

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
if (
    st.session_state.analyzed
    and option not in [
        "Upload Multiple FASTA",
        "Multiple Sequence Alignment",
        "Pairwise Sequence Alignment"
    ]
):
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