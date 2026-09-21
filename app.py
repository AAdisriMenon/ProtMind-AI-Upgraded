import streamlit as st
from backend import apply_custom_css

st.set_page_config(page_title="ProtMind AI Dashboard", page_icon="🧠", layout="centered")
apply_custom_css()

# Initialize Session State Memory to pass data between pages
if 'payload' not in st.session_state:
    st.session_state['payload'] = None

st.markdown("<h1>🧠 ProtMind AI</h1>", unsafe_allow_html=True)
st.markdown("<h3>Explainable Protein Mutation Intelligence & Therapeutic Discovery</h3>", unsafe_allow_html=True)
st.divider()

st.info("👋 Welcome to the ProtMind AI Pipeline. Please use the sidebar to navigate through the modules, starting with Step 1: User Input.")

if st.session_state['payload']:
    st.success(f"✅ Active Session Found: Analyzing {st.session_state['payload']['uniprot_id']} ({st.session_state['payload']['mutation']})")

    results = st.session_state.get('results', {})
    stage_map = [
        ("Preprocessing (3)", "features"), ("Mutation Analysis (4)", "pathogenicity"),
        ("Stability (6)", "stability"), ("Stereochemistry (7)", "angles"),
        ("Explainable AI (8)", "composite"), ("Therapeutic Targeting (9)", "drugs"),
    ]
    done = sum(1 for _, key in stage_map if key in results)
    st.progress(done / len(stage_map), text=f"Pipeline progress: {done}/{len(stage_map)} analysis stages completed this session")
st.markdown("### Stereochemical Validation")
st.markdown("Evaluating backbone strain and torsional feasibility of the mutated residue.")

# Render your plot
fig = generate_ramachandran_plot(phi, psi, mutation)
st.plotly_chart(fig, use_container_width=True)

# Description Box matching your layout design
st.info("""
**Angles falling in the blank (white) areas of the plot represent steric clashes—the atoms are physically overlapping.**
""")

st.markdown("""
- **Allowed (Pink):** Optimal backbone torsion without steric clashes.
- **Partially Allowed (Blue):** Torsion permitted with slight conformational strain.
- **Disallowed (White):** Severe steric overlap; mutation is highly unstable.
""")
