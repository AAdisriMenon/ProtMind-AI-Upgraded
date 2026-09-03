import streamlit as st
import time
import py3Dmol
from stmol import showmol
from backend import apply_custom_css, fetch_alphafold_url, fetch_pdb_content

st.set_page_config(page_title="Step 5 | ProtMind AI", page_icon="🧩")
apply_custom_css()

st.markdown("<h2>🧩 Step 5: Protein Structure Prediction</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']
uniprot_id = payload['uniprot_id']
mutation = payload['mutation']

try:
    mut_pos = int(mutation[1:-1])
except ValueError:
    mut_pos = None

with st.spinner(f"🔬 Fetching AlphaFold 3D structure for {uniprot_id}..."):
    af_pdb_url = fetch_alphafold_url(uniprot_id)
    pdb_content = fetch_pdb_content(af_pdb_url)
    time.sleep(1)

if not pdb_content:
    st.error(f"❌ Could not retrieve 3D structure for {uniprot_id}. It may not be available in the AlphaFold database.")
else:
    col1, col2 = st.columns([2, 1])

    with col2:
        st.markdown("<h3>📊 Visualization Controls</h3>", unsafe_allow_html=True)
        
        # Interactive dropdowns for structural views
        view_style = st.selectbox(
            "🎨 Select Protein Style", 
            ["Cartoon (Ribbon)", "Stick", "Sphere (Space-filling)", "Line"]
        )
        
        color_scheme = st.selectbox(
            "🌈 Select Color Scheme",
            ["Default (Cyan)", "By Secondary Structure", "Spectrum (Rainbow)"]
        )
        
        st.markdown("---")
        st.info("The target mutated residue is always highlighted in **RED** so you don't lose it when changing styles.")
        
        st.metric(label="Structure Source", value="AlphaFold2")
        st.metric(label="Target Mutation", value=mutation)
        
        st.markdown("<br>", unsafe_allow_html=True)
        st.download_button(
            label="📥 Download .PDB File",
            data=pdb_content,
            file_name=f"{uniprot_id}_AF.pdb",
            mime="text/plain"
        )

    with col1:
        st.markdown("<h3>🔍 Interactive 3D Viewer</h3>", unsafe_allow_html=True)
        
        # Initialize py3Dmol Viewer
        view = py3Dmol.view(width=500, height=500)
        view.addModel(pdb_content, 'pdb')
        view.setBackgroundColor('#161A23')
        
        # 1. Determine the color configuration based on user selection
        if color_scheme == "By Secondary Structure":
            color_config = {'colorscheme': 'ssPyMOL'}
        elif color_scheme == "Spectrum (Rainbow)":
            color_config = {'colorscheme': 'amino'}
        else:
            color_config = {'color': '#00F5D4'}
            
        # 2. Apply the base structural style based on user selection
        if view_style == "Cartoon (Ribbon)":
            view.setStyle({'cartoon': color_config})
        elif view_style == "Stick":
            view.setStyle({'stick': color_config})
        elif view_style == "Sphere (Space-filling)":
            view.setStyle({'sphere': color_config})
        elif view_style == "Line":
            view.setStyle({'line': color_config})
        
        # 3. Always highlight the specific mutated residue in Red on top of the base style
        if mut_pos:
            view.addStyle({'resi': str(mut_pos)}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.25}})
            view.addStyle({'resi': str(mut_pos)}, {'sphere': {'color': '#FF0055', 'radius': 1.2}})
            view.addLabel(
                f"Mut: {mutation}", 
                {'fontOpacity': 1, 'fontSize': 14, 'fontColor': 'white', 'backgroundColor': '#FF0055'}, 
                {'resi': str(mut_pos)}
            )
            # Auto-zoom directly to the mutation site
            view.zoomTo({'resi': str(mut_pos)})
        else:
            view.zoomTo()
            
        # Render the py3Dmol object in Streamlit
        showmol(view, height=500, width=500)