import streamlit as st
from backend import apply_custom_css, validate_sequence, validate_uniprot_id, validate_mutation, parse_fasta

st.set_page_config(page_title="Step 1 | ProtMind AI", page_icon="⚡")
apply_custom_css()

st.markdown("<h2>⚡ Step 1: User Input</h2>", unsafe_allow_html=True)

with st.form("user_input_form", clear_on_submit=False):
    col1, col2 = st.columns(2)
    with col1:
        uniprot_id = st.text_input("🧬 UniProt ID", placeholder="P04637", help="Example: P04637")
    with col2:
        mutation = st.text_input("🧪 Mutation", placeholder="R175H", help="Format: [WT][Position][Mut] (e.g., R175H)")

    st.markdown("---")
    st.markdown("### 📂 Sequence Input Options")

    # Option A: File Upload
    uploaded_file = st.file_uploader(
        "Upload FASTA File (.fasta, .fa, .txt)",
        type=["fasta", "fa", "txt"],
        help="Upload a FASTA file containing your target protein sequence."
    )

    st.markdown("<p style='text-align:center;'>— OR —</p>", unsafe_allow_html=True)

    # Option B: Manual Sequence Paste
    sequence_input = st.text_area(
        "Paste Protein Sequence or FASTA text",
        placeholder=">sp|P04637|...\nMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPL...",
        height=150
    )

    submitted = st.form_submit_button("🚀 Validate & Save Data")

if submitted:
    errors = []
    final_sequence = ""
    source_type = ""

    # Priority 1: File Upload
    if uploaded_file is not None:
        try:
            raw_text = uploaded_file.getvalue().decode("utf-8")
            final_sequence = parse_fasta(raw_text)
            source_type = f"Uploaded file ({uploaded_file.name})"
        except Exception:
            errors.append("❌ Could not read the uploaded FASTA file. Ensure it is UTF-8 encoded text.")
    # Priority 2: Text Area
    elif sequence_input.strip():
        final_sequence = parse_fasta(sequence_input)
        source_type = "Manual input"
    else:
        errors.append("❌ Please upload a FASTA file or paste a protein sequence.")

    # Validation Checks
    if final_sequence and not validate_sequence(final_sequence):
        errors.append("❌ Invalid amino acid characters detected. Only standard 20 amino acids allowed.")

    if not uniprot_id or not validate_uniprot_id(uniprot_id):
        errors.append("❌ Invalid UniProt ID format. Example: P04637")

    if not mutation or not validate_mutation(mutation):
        errors.append("❌ Invalid Mutation format. Example: R175H")

    # Output Handling
    if errors:
        for error in errors:
            st.error(error)
    else:
        st.session_state['payload'] = {
            "sequence": final_sequence.strip().upper(),
            "uniprot_id": uniprot_id.strip().upper(),
            "mutation": mutation.strip().upper(),
            "source": source_type
        }
        st.success(f"✅ Sequence successfully extracted and validated via {source_type}!")
        st.info(f"**Loaded Sequence Length:** {len(final_sequence)} amino acids")
        st.write("👉 Proceed to **Step 2 (Data Retrieval)** using the sidebar navigation.")