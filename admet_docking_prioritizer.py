# The ADMET & Docking Prioritizer App (Streamlit/RDKit Implementation)
# ALL LOGIC HAS BEEN MOVED INTO THE MAIN EXECUTION FLOW (WITHOUT 'def')

import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from io import StringIO
# Note: Draw is imported but not used in the final logic

# --- 1. Streamlit UI Setup ---

# Set Streamlit page configuration
st.set_page_config(layout="wide", page_title="CADD Prioritizer")

st.title("🔬 ADMET & Docking Prioritizer App")
st.markdown("A tool demonstrating the synergy between Medicinal Chemistry, Cheminformatics, and Docking for compound prioritization.")

# Example Data for User Input
EXAMPLE_DATA = """SMILES,Docking_Score
CC(=O)Oc1ccccc1C(=O)O,-7.2 # Aspirin (Drug-like, good score)
COc1ccc(C(C)Nc2ccc(C)cc2)cc1,-6.5 # Example Pass
O=C1CCc2c(C)nc(C)c2N1C1CC1,-4.1 # Example Fail (High LogP)
CC(C)CN(C)CC(O)C(C)c1ccc(O)c(Cl)c1,-5.9 # Example Fail (Many H Acceptors)
C1=CC=C2C(=C1)C=CC(=O)C2=O,-7.5 # Example Pass (Very small, good score)
InvalidSMILES,-8.0 # Intentional error
"""

# User Input Area
st.header("Input Candidate Molecules")
input_method = st.radio(
    "Select Input Method:",
    ('Use Example Data', 'Paste Custom Data (SMILES, Docking_Score)'),
    horizontal=True
)

# Assign raw data based on selection
if input_method == 'Use Example Data':
    raw_data = EXAMPLE_DATA
else:
    raw_data = st.text_area(
        "Paste Data (CSV format: SMILES,Docking_Score)",
        value=EXAMPLE_DATA,
        height=250
    )

# Process Button
if st.button("Analyze & Prioritize Candidates", type="primary"):
    
    try:
        # Load data into DataFrame
        df_input = pd.read_csv(StringIO(raw_data.strip()))
        
        if 'SMILES' not in df_input.columns or 'Docking_Score' not in df_input.columns:
            st.error("Input data must contain 'SMILES' and 'Docking_Score' columns.")
        else:
            with st.spinner('Calculating properties and applying filters...'):
                
                # --- CORE LOGIC: No functions used ---
                results = []
                
                for index, row in df_input.iterrows():
                    smiles = row['SMILES']
                    docking_score = row['Docking_Score']

                    # Programming & Cheminformatics: RDKit Conversion
                    mol = Chem.MolFromSmiles(smiles)

                    if mol is None:
                        # Handle invalid SMILES input
                        results.append({
                            'SMILES': smiles,
                            'Docking_Score': docking_score,
                            'Status': 'Invalid SMILES',
                            'MW': None, 'LogP': None, 'HDonors': None, 'HAcceptors': None, 'Violations': None
                        })
                        continue

                    # Foundational Sciences: Drug-Likeness Check & Property Calculation
                    
                    # Calculate RDKit Descriptors
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    h_donors = Descriptors.NumHDonors(mol)
                    h_acceptors = Descriptors.NumHAcceptors(mol)
                    
                    # Apply Lipinski's Rule of Five Logic
                    violations = 0
                    if mw > 500: violations += 1
                    if logp > 5: violations += 1
                    if h_donors > 5: violations += 1
                    if h_acceptors > 10: violations += 1
                    
                    is_drug_like = (violations <= 1)
                    
                    status = 'Pass' if is_drug_like else 'Fail (Lipinski Violation)'
                    
                    results.append({
                        'SMILES': smiles,
                        'Docking_Score': docking_score, # Specialized Technique metric
                        'Status': status,
                        'MW': round(mw, 2),
                        'LogP': round(logp, 2),
                        'HDonors': h_donors,
                        'HAcceptors': h_acceptors,
                        'Violations': violations
                    })

                # Convert results to DataFrame
                df_results = pd.DataFrame(results)
                
                # Affinity Prioritization (Specialized Technique): Rank candidates that passed the filter
                
                # Separate passing and failing molecules
                df_pass = df_results[df_results['Status'] == 'Pass'].copy()
                df_fail = df_results[df_results['Status'] != 'Pass'].copy()

                # Rank 'Pass' molecules based on Docking_Score (ascending for lowest score = rank 1)
                df_pass.sort_values(by='Docking_Score', ascending=True, inplace=True)
                df_pass['Final_Rank'] = range(1, len(df_pass) + 1)
                
                # 'Fail' molecules are not ranked
                df_fail['Final_Rank'] = '-'
                
                # Recombine and sort for final output
                df_final = pd.concat([df_pass, df_fail]).sort_values(by=['Status', 'Final_Rank'], 
                                                                  ascending=[False, True], 
                                                                  ignore_index=True)
                
                # --- END OF CORE LOGIC ---
                
            st.success("Analysis Complete!")
            
            st.header("Results: Prioritized Drug Candidates")
            st.caption("Lower Docking Scores indicate higher predicted affinity. Only 'Pass' candidates are ranked.")
            
            # --- Visualization ---

            # Note: Dynamic coloring via applymap was removed to avoid using the 'def' keyword.
            st.dataframe(
                df_final,
                use_container_width=True,
                column_config={
                    "Docking_Score": st.column_config.NumberColumn("Docking Score (kcal/mol)", format="%.2f"),
                    "MW": st.column_config.NumberColumn("MW", format="%.2f"),
                    "LogP": st.column_config.NumberColumn("LogP", format="%.2f"),
                    "HDonors": st.column_config.NumberColumn("H Donors"),
                    "HAcceptors": st.column_config.NumberColumn("H Acceptors"),
                    "Violations": st.column_config.NumberColumn("Rof5 Violations"),
                    "Final_Rank": st.column_config.TextColumn("Final Rank"),
                }
            )
            
            # Additional Section: Drug-Likeness Summary
            pass_count = (df_final['Status'] == 'Pass').sum()
            fail_count = (df_final['Status'].str.contains('Fail')).sum()
            st.subheader("Summary")
            st.info(f"**{pass_count}** candidate(s) passed the Drug-Likeness Filter and were prioritized by Docking Score.")
            st.warning(f"**{fail_count}** candidate(s) failed the Drug-Likeness Filter and were excluded from ranking.")
            
            st.markdown("---")
            st.markdown("### Technical Breakdown (The Prerequisites in Action)")
            st.code(
                f"""
                Total Candidates: {len(df_input)}
                
                # Foundational Sciences (Medicinal Chemistry) & RDKit (Cheminformatics):
                # Lipinski's Rule of Five applied:
                # MW <= 500
                # LogP <= 5
                # H Donors <= 5
                # H Acceptors <= 10

                # Specialized Techniques (Docking):
                # Priority ranking based on 'Docking_Score' (lowest score is best) 
                # only applied to candidates that passed the Medicinal Chemistry filter.
                """
            )

    except Exception as e:
        st.error(f"An error occurred during data processing: {e}")

# Initial Instructions if no button is clicked
else:
    st.info("Click 'Analyze & Prioritize Candidates' to run the CADD workflow on the input data.")

