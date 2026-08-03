import streamlit as st
import plotly.express as px
import pandas as pd

from modules.analysis import (
    sequence_length,
    nucleotide_count,
    gc_content,
    at_content,
    reverse_complement,
    transcribe,
    translate,
    find_orf,
    motif_search,
    restriction_sites,
    amino_acid_frequency,
    colored_sequence,
    codon_frequency
)
def analyze_sequence(dna, motif, unique_key=""):

    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs([
        "📊 Basic",
        "🧬 RNA",
        "🧬 ORF",
        "🔍 Motif",
        "✂ Restriction",
        "🧪 Amino",
        "🎨 DNA",
        "🧬 Codons",
        "📄 Report"
    ])

    # ==========================
    # TAB 1
    # ==========================

    with tab1:

        st.subheader("📊 Basic Analysis")

        col1, col2, col3 = st.columns(3)

        col1.metric("Length", sequence_length(dna))
        col2.metric("GC %", gc_content(dna))
        col3.metric("AT %", at_content(dna))

        counts = nucleotide_count(dna)

        df_counts = pd.DataFrame({
            "Nucleotide": ["A", "T", "G", "C"],
            "Count": [
                counts["A"],
                counts["T"],
                counts["G"],
                counts["C"]
            ]
        })

        st.dataframe(df_counts, width="stretch")

        fig = px.bar(
            df_counts,
            x="Nucleotide",
            y="Count",
            text="Count",
            title="Nucleotide Distribution"
        )

        st.plotly_chart(
            fig,
            width="stretch",
            key=f"basic_{unique_key}_{len(dna)}"
        )

        pie_df = pd.DataFrame({
            "Content": ["GC", "AT"],
            "Percentage": [
                gc_content(dna),
                at_content(dna)
            ]
        })

        pie = px.pie(
            pie_df,
            names="Content",
            values="Percentage",
            title="GC vs AT Composition"
        )

        st.plotly_chart(
            pie,
            width="stretch",
            key=f"pie_{unique_key}_{len(dna)}"
        )

    # ==========================
    # TAB 2
    # ==========================

    with tab2:

        st.subheader("Reverse Complement")
        st.code(reverse_complement(dna))

        st.subheader("RNA Transcription")
        st.code(transcribe(dna))

        st.subheader("Protein Translation")
        st.code(translate(dna))
    # ==========================
    # TAB 3
    # ==========================

    with tab3:

        st.subheader("🧬 Open Reading Frame (ORF)")

        st.code(find_orf(dna))

    # ==========================
    # TAB 4
    # ==========================

    with tab4:

        st.subheader("🔍 Motif Search")

        if motif:

            positions = motif_search(dna, motif)

            if positions:
                st.success(f"Motif Found at: {positions}")
            else:
                st.warning("Motif Not Found")

        else:
            st.info("Enter a motif above.")

    # ==========================
    # TAB 5
    # ==========================

    with tab5:

        st.subheader("✂ Restriction Enzyme Analysis")

        sites = restriction_sites(dna)

        for enzyme, positions in sites.items():

            if positions:
                st.success(f"{enzyme}: {positions}")
            else:
                st.info(f"{enzyme}: Not Found")

    # ==========================
    # TAB 6
    # ==========================

    with tab6:

        st.subheader("🧪 Amino Acid Frequency")

        aa = amino_acid_frequency(dna)

        if aa:

            aa_df = pd.DataFrame({
                "Amino Acid": list(aa.keys()),
                "Count": list(aa.values())
            })

            st.dataframe(aa_df, width="stretch")

            fig = px.bar(
                aa_df,
                x="Amino Acid",
                y="Count",
                text="Count",
                title="Amino Acid Frequency"
            )

            st.plotly_chart(
                fig,
                width="stretch",
                key=f"aa_{unique_key}_{len(dna)}"
            )

        else:

            st.info("No amino acids found.")

    # ==========================
    # TAB 7
    # ==========================

    with tab7:

        st.subheader("🎨 Colored DNA Sequence")

        if len(dna) > 1000:
            st.warning("Sequence is very large. Showing first 1000 nucleotides only.")
            display_seq = dna[:1000]
        else:
            display_seq = dna

        st.markdown(
            colored_sequence(display_seq),
            unsafe_allow_html=True
        )

    # ==========================
    # TAB 8
    # ==========================

    with tab8:

        st.subheader("🧬 Codon Frequency")

        codons = codon_frequency(dna)

        if codons:

            codon_df = pd.DataFrame({
                "Codon": list(codons.keys()),
                "Count": list(codons.values())
            })

            st.dataframe(codon_df, width="stretch")

            fig = px.bar(
                codon_df,
                x="Codon",
                y="Count",
                text="Count",
                title="Codon Frequency"
            )

            st.plotly_chart(
                fig,
                width="stretch",
                key=f"codon_{unique_key}_{len(dna)}"
            )

        else:

            st.info("No codons found.")

    # ==========================
    # REPORT
    # ==========================

    report = f"""
DNA Sequence Analysis Report
============================

Sequence:
{dna[:1000]}

Length: {sequence_length(dna)}

GC Content: {gc_content(dna)}%

AT Content: {at_content(dna)}%

Nucleotide Counts
A : {counts['A']}
T : {counts['T']}
G : {counts['G']}
C : {counts['C']}

Reverse Complement
{reverse_complement(dna)}

RNA
{transcribe(dna)}

Protein
{translate(dna)}

ORF
{find_orf(dna)}
"""

    # ==========================
    # TAB 9
    # ==========================

    with tab9:

        st.subheader("📄 Download Report")

        st.download_button(
            label="📥 Download Report",
            data=report,
            file_name="DNA_Report.txt",
            mime="text/plain",
            key=f"report_{unique_key}_{len(dna)}"
        )