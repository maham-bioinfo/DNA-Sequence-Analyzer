import streamlit as st
import plotly.express as px
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go


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
    codon_frequency,
    gc_sliding_window,
    dna_quality_score,
    genome_features,
    genome_map_data,
    genome_dataframe,
    pairwise_alignment
)
def analyze_sequence(dna, motif, unique_key=""):

    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10, tab11, tab12 = st.tabs([
        "📊 Basic",
        "🧬 RNA",
        "🧬 ORF",
        "🔍 Motif",
        "✂ Restriction",
        "🧪 Amino",
        "🎨 DNA",
        "🧬 Codons",
        "📄 Report",
        "📈 GC Graph",
        "⭐ Quality",
        "🗺 Genome Browser"
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
    # ==========================
    # TAB 10
    # ==========================

    with tab10:

        st.subheader("📈 GC Content Sliding Window")

        window = st.slider(
            "Window Size",
            20,
            500,
            100,
            key=f"window_{unique_key}"
        )

        positions, gc = gc_sliding_window(
            dna,
            window
        )

        gc_df = pd.DataFrame({
            "Position": positions,
            "GC %": gc
        })

        fig = px.line(
            gc_df,
            x="Position",
            y="GC %",
            title="GC Content Across Genome"
        )

        st.plotly_chart(
            fig,
            width="stretch",
            key=f"gc_graph_{unique_key}"
        )
    # ==========================
    # TAB 11
    # ==========================

    with tab11:

        st.subheader("⭐ DNA Quality Score")

        score, report = dna_quality_score(dna)

        st.metric(
            "Overall Score",
            f"{score}/100"
        )

        progress = score / 100

        st.progress(progress)

        st.divider()

        for text, status in report:

            if status:
                st.success(text)
            else:
                st.warning(text)

    # ==========================
    # TAB 12
    # ==========================

    with tab12:

        st.subheader("🗺 Interactive Genome Browser")

        df = genome_dataframe(dna, motif)

        if df.empty:

            st.info("No genomic features found.")

        else:

            show_labels = st.toggle(
                "Show Labels",
                value=True,
                key=f"show_labels_{unique_key}"
            )

            fig = go.Figure()

            color_map = {

                "ORF":"green",

                "Motif":"orange",

                "EcoRI":"red",

                "BamHI":"blue",

                "HindIII":"purple",

                "NotI":"brown"

            }

            marker_map = {

                "ORF":"square",

                "Motif":"diamond",

                "EcoRI":"triangle-up",

                "BamHI":"triangle-down",

                "HindIII":"circle",

                "NotI":"x"

            }
            for _, row in df.iterrows():

                feature = row["Feature"]

                color = color_map.get(feature, "black")

                symbol = marker_map.get(feature, "circle")

                if feature == "ORF":

                    fig.add_shape(

                        type="rect",

                        x0=row["Start"],

                        x1=row["End"],

                        y0=row["Track"]-0.25,

                        y1=row["Track"]+0.25,

                        line=dict(color=color),

                        fillcolor=color,

                        opacity=0.75

                    )

                    if show_labels:

                        fig.add_annotation(

                            x=(row["Start"]+row["End"])/2,

                            y=row["Track"]+0.45,

                            text="ORF",

                            showarrow=False,

                            font=dict(size=10,color="green")

                        )

                else:

                    fig.add_trace(

                        go.Scatter(

                            x=[row["Start"]],

                            y=[row["Track"]],

                            mode="markers+text" if show_labels else "markers",

                            text=[feature] if show_labels else None,

                            textposition="top center",

                            marker=dict(

                                size=12,

                                color=color,

                                symbol=symbol

                            ),

                            hovertemplate=(

                                f"<b>{feature}</b><br>"

                                f"Start : {row['Start']}<br>"

                                f"End : {row['End']}<extra></extra>"

                            ),

                            showlegend=False

                        )

                    )
            fig.update_layout(

                title="Interactive Genome Browser",

                template="plotly_white",

                height=550,

                hovermode="closest",

                dragmode="pan",

                plot_bgcolor="white",

                paper_bgcolor="white",

                margin=dict(

                    l=70,

                    r=40,

                    t=70,

                    b=60

                )

            )

            fig.update_xaxes(

                title="Genome Position (bp)",

                showgrid=True,

                rangeslider=dict(

                    visible=True

                ),

                zeroline=False

            )

            fig.update_yaxes(

                title="Genome Features",

                tickmode="array",

                tickvals=[1,2,3,4,5,6],

                ticktext=[

                    "ORF",

                    "Motif",

                    "EcoRI",

                    "BamHI",

                    "HindIII",

                    "NotI"

                ],

                showgrid=True,

                autorange="reversed"

            )

            st.plotly_chart(

                fig,

                width="stretch",
                key=f"genome_browser_{unique_key}",
                config={

                    "displaylogo": False,

                    "scrollZoom": True,

                    "modeBarButtonsToAdd":[

                        "zoom2d",

                        "pan2d",

                        "select2d",

                        "lasso2d",

                        "zoomIn2d",

                        "zoomOut2d",

                        "resetScale2d"

                    ]

                }

            )
def show_pairwise_alignment(alignment):

    st.subheader("🧬 Pairwise Sequence Alignment")

    st.metric(
        "Alignment Score",
        round(alignment.score,2)
    )

    st.code(str(alignment))
def show_msa_alignment(alignment):

    st.subheader(
        "🧬 Multiple Sequence Alignment"
    )

    # ---------------------------------
    # Alignment statistics
    # ---------------------------------

    num_sequences = len(alignment)

    alignment_length = alignment.get_alignment_length()

    col1, col2 = st.columns(2)

    with col1:

        st.metric(
            "Sequences",
            num_sequences
        )

    with col2:

        st.metric(
            "Alignment Length",
            f"{alignment_length} bp"
        )

    st.divider()

    # ---------------------------------
    # Display alignment
    # ---------------------------------

    st.subheader(
        "Aligned Sequences"
    )

    for record in alignment:

        st.markdown(
            f"**{record.id}**"
        )

        st.code(
            str(record.seq),
            language="text"
        )

    # ---------------------------------
    # Download aligned FASTA
    # ---------------------------------

    from io import StringIO
    from Bio import AlignIO

    output = StringIO()

    AlignIO.write(
        alignment,
        output,
        "fasta"
    )

    aligned_fasta = output.getvalue()

    st.download_button(
        "📥 Download Aligned FASTA",
        aligned_fasta,
        file_name="aligned_sequences.fasta",
        mime="text/plain"
    )