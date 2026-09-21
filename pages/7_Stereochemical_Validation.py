import streamlit as st
import time
from backend import apply_custom_css, fetch_alphafold_url, fetch_pdb_content, calculate_phi_psi, generate_ramachandran_plot

st.set_page_config(page_title="Step 7 | ProtMind AI", page_icon="📐")
apply_custom_css()

st.markdown("<h2>📐 Step 7: Stereochemical Validation</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']
mutation = payload['mutation']
mut_pos = int(mutation[1:-1])

with st.spinner("📐 Parsing AlphaFold atomic coordinates for Phi/Psi angles..."):
    af_pdb_url = fetch_alphafold_url(payload['uniprot_id'])
    pdb_content = fetch_pdb_content(af_pdb_url)
    angles = calculate_phi_psi(pdb_content, mut_pos)
    time.sleep(1)

st.session_state.setdefault('results', {})
st.session_state['results']['angles'] = angles

if not angles:
    st.error("❌ Could not calculate backbone angles. The residue position might not exist in the AlphaFold structure.")
else:
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.markdown("<h3>Backbone Torsion</h3>", unsafe_allow_html=True)
        st.metric(label="Phi (Φ) Angle", value=f"{angles['phi']:.2f}°")
        st.metric(label="Psi (Ψ) Angle", value=f"{angles['psi']:.2f}°")
        st.info("Angles falling in the blank (dark) areas of the plot represent steric clashes—the atoms are physically overlapping.")
        st.info("Angles falling in the blank (white) areas of the plot represent steric clashes—the atoms are physically overlapping.")
    
        st.markdown("""
        - **Allowed (Pink):** Optimal backbone torsion without steric clashes.
        - **Partially Allowed (Blue):** Torsion permitted with slight conformational strain.
        - **Disallowed (White):** Severe steric overlap; mutation is highly unstable.
        """)

        
    with col2:
        fig = generate_ramachandran_plot(angles['phi'], angles['psi'], mutation)
        st.plotly_chart(fig, use_container_width=True)
