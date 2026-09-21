import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def classify_ramachandran(phi, psi):
    """Classifies the structural region and tolerance of the dihedral angles."""
    # Top-Left: Beta-sheet
    if -180 <= phi <= -20 and 30 <= psi <= 180:
        if -150 <= phi <= -50 and 90 <= psi <= 180:
            return "Beta-sheet", "✅ Allowed (Optimal Sterics)"
        return "Beta-sheet", "⚠️ Partially Allowed"
        
    # Bottom-Left: Right-handed alpha-helix
    elif -180 <= phi <= -20 and -180 <= psi <= 30:
        if -140 <= phi <= -50 and -90 <= psi <= -30:
            return "Right-handed α-helix", "✅ Allowed (Optimal Sterics)"
        return "Right-handed α-helix", "⚠️ Partially Allowed"
        
    # Top-Right: Left-handed alpha-helix
    elif 20 <= phi <= 100 and 20 <= psi <= 100:
        return "Left-handed α-helix", "⚠️ Partially Allowed"
        
    # Everywhere else
    else:
        return "Unstructured / Loop", "❌ Disallowed (High Steric Clash)"

def plot_ramachandran_reference(phi, psi, label="Mutant"):
    """Recreates the specific Ramachandran plot style from your reference image."""
    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Base styling matching the reference (White background, green text)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Draw crosshairs
    ax.axhline(0, color='black', linewidth=1.5)
    ax.axvline(0, color='black', linewidth=1.5)
    
    # --- Define Polygons for the Regions ---
    # Colors matching reference image
    c_partial = '#63d3ff' # Light Blue
    c_allowed = '#cc338b' # Pink/Magenta
    
    # 1. Beta-Sheet Regions (Top Left)
    beta_blue = patches.Polygon([(-180, 20), (-180, 180), (-40, 180), (-20, 130), (-40, 80), (-80, 20)], facecolor=c_partial, edgecolor='black', linewidth=1.2)
    beta_pink = patches.Polygon([(-150, 90), (-150, 180), (-60, 180), (-40, 140), (-60, 90)], facecolor=c_allowed, edgecolor='black', linewidth=1.2)
    
    # 2. Right-Handed Helix Regions (Bottom Left)
    rh_blue = patches.Polygon([(-180, -180), (-40, -180), (-20, -130), (-20, 0), (-50, 20), (-180, 20)], facecolor=c_partial, edgecolor='black', linewidth=1.2)
    rh_pink = patches.Polygon([(-150, -180), (-50, -180), (-40, -140), (-40, -30), (-70, -30), (-100, -60), (-150, -60)], facecolor=c_allowed, edgecolor='black', linewidth=1.2)
    
    # 3. Left-Handed Helix Region (Right Side)
    lh_blue = patches.Polygon([(30, 20), (30, 100), (80, 100), (80, 20)], facecolor=c_partial, edgecolor='black', linewidth=1.2)

    # Add patches to plot
    ax.add_patch(beta_blue)
    ax.add_patch(rh_blue)
    ax.add_patch(lh_blue)
    ax.add_patch(beta_pink)
    ax.add_patch(rh_pink)
    
    # Add text labels inside regions matching reference
    ax.text(-110, 130, 'Beta-sheet', color='white', weight='bold', fontsize=10, ha='center')
    ax.text(-100, -60, 'Right handed\nalpha-helix.', color='#008000', weight='bold', fontsize=9, ha='center', va='top')
    ax.text(90, 60, 'Left\nhanded\nalpha-helix.', color='#008000', weight='bold', fontsize=9, ha='left', va='center')
    
    # Plot the specific mutation angle
    ax.scatter(phi, psi, color='yellow', s=180, edgecolor='black', linewidth=2, zorder=5, label=f"{label} ($\phi$={phi:.1f}, $\psi$={psi:.1f})")
    
    # Axis formatting to match reference image exactly
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks([-180, 0, 180])
    ax.set_yticks([-180, 0, 180])
    
    ax.set_xticklabels(['-180', '0', '+180'], color='#00aa00', weight='bold', fontsize=12)
    ax.set_yticklabels(['-180', '0', '+180'], color='#00aa00', weight='bold', fontsize=12)
    
    # Custom axis labels placed like the reference image
    ax.text(-90, -200, '-phi', color='#00aa00', weight='bold', fontsize=12, ha='center')
    ax.text(90, -200, '+phi', color='#00aa00', weight='bold', fontsize=12, ha='center')
    ax.text(-210, -90, '-psi', color='#00aa00', weight='bold', fontsize=12, va='center')
    ax.text(-210, 90, '+psi', color='#00aa00', weight='bold', fontsize=12, va='center')
    
    ax.set_title("The Ramachandran Plot", color='#00aa00', weight='bold', fontsize=14, pad=15)
    
    # Remove standard spines, use only outer box
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        
    ax.legend(loc='upper right', framealpha=1, facecolor='white', edgecolor='black')
    
    return fig

# ==========================================
# STREAMLIT UI INTEGRATION
# ==========================================

st.markdown("### 🧬 Stereochemical Validation")
st.write("Evaluating backbone strain and torsional feasibility of the mutated residue.")

# Fetch angles from your session state
angles = st.session_state.get('results', {}).get('angles', None)
payload = st.session_state.get('payload', {})

if angles and 'phi' in angles and 'psi' in angles:
    phi_val = angles['phi']
    psi_val = angles['psi']
    
    wt = payload.get('mutation', 'X')[0]
    pos = payload.get('mutation', '0X')[1:-1]
    mut = payload.get('mutation', 'X')[-1]
    
    # Get Classification
    structure_type, sterics_status = classify_ramachandran(phi_val, psi_val)
    
    col1, col2 = st.columns([1.5, 1])
    
    with col1:
        # Generate and display the plot
        fig = plot_ramachandran_reference(phi_val, psi_val, label=f"{mut}{pos}")
        st.pyplot(fig)
        
    with col2:
        st.markdown("#### Structural Diagnostics")
        st.info(f"**$\phi$ (Phi) Angle:** {phi_val:.2f}°")
        st.info(f"**$\psi$ (Psi) Angle:** {psi_val:.2f}°")
        
        st.markdown("#### Classification")
        st.success(f"**Secondary Structure:** {structure_type}")
        
        if "Allowed" in sterics_status and "Partially" not in sterics_status:
            st.success(f"**Steric Feasibility:** {sterics_status}")
        elif "Partially" in sterics_status:
            st.warning(f"**Steric Feasibility:** {sterics_status}")
        else:
            st.error(f"**Steric Feasibility:** {sterics_status}")
            
        st.markdown("""
        * **Allowed (Pink):** Optimal backbone torsion without steric clashes.
        * **Partially Allowed (Blue):** Torsion permitted with slight conformational strain.
        * **Disallowed (White):** Severe steric overlap; mutation is highly unstable.
        """)
else:
    st.warning("Dihedral angles ($\phi$, $\psi$) not found. Ensure 3D structures were processed in Step 4.")
