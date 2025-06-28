# for_fragpipe.py - Modified version with function interface

import re
import operator
import numpy as np
import pandas as pd
import os
import sys
import irt_alignment as irt_align

def merge_libraries(sample_library_path, systemhc_lib_paths, output_dir):
    """
    Merge sample library with SysteMHC libraries
    
    Parameters:
    -----------
    sample_library_path : str
        Path to the sample library TSV file
    systemhc_lib_paths : list
        List of paths to SysteMHC library TSV files
    output_dir : str
        Directory where the merged library will be saved
    
    Returns:
    --------
    str
        Path to the merged library file
    """
    print(f"Sample library: {sample_library_path}")
    print(f"SysteMHC libraries: {systemhc_lib_paths}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load sample library
    sample_library = pd.read_csv(sample_library_path, sep='\t')
    sample_library2 = sample_library.copy()
    sample_library2['ions'] = sample_library2['ModifiedPeptideSequence'] + sample_library2['PrecursorCharge'].astype(str)
    
    # Load irt data for RT normalization (if available)
    try:
        # Try to use the irt file in the output directory
        irt_file_path = os.path.join(output_dir, 'irt_SYSTEMHC.csv')
        if os.path.exists(irt_file_path):
            df_need = pd.read_csv(irt_file_path)
        else:
            # if file doesn't exist in output dir
            print(f"Failed:  {irt_file_path} not found")
        
        # RT normalization
        rt_reference_run = sample_library[['ModifiedPeptideSequence','PrecursorCharge','NormalizedRetentionTime']].drop_duplicates()
        rt_reference_run.columns = ['modified_peptide','precursor_charge','irt']
    
        df_need.columns = ['modified_peptide','precursor_charge','RT']
    
        aligned_runs1 = irt_align.lowess2(df_need, rt_reference_run, 'RT', 'irt', 0.01, 0, 10)
        pepida1 = aligned_runs1
        pepida1 = pepida1.loc[np.isfinite(pepida1['irt'])]
        pqp = pepida1
        pqp2 = pqp.groupby(['modified_peptide','precursor_charge'])[['irt']].median().reset_index()
        pqp2.columns = ['ModifiedPeptide','PrecursorCharge','NormalizedRetentionTime']
    
        # Save RT alignment results
        rt_aligned_path = os.path.join(output_dir, 'rt_aligned2reference.csv')
        pqp2.to_csv(rt_aligned_path, index=False)
        print(f"RT alignment results saved to: {rt_aligned_path}")
    except Exception as e:
        print(f"Warning: RT normalization skipped - {str(e)}")
        # Continue without RT normalization
    
    # Combine SysteMHC libraries
    systemhc_libs = []
    for libp in systemhc_lib_paths:
        try:
            datmp = pd.read_csv(libp, sep='\t')
            systemhc_libs.append(datmp)
            print(f"Loaded SysteMHC library: {libp}")
        except Exception as e:
            print(f"Warning: Failed to load {libp} - {str(e)}")
    
    if not systemhc_libs:
        raise Exception("No valid SysteMHC libraries were loaded")
    
    da = pd.concat(systemhc_libs)
    da = da.drop_duplicates()
    da1 = da.copy()
    da1['ions'] = da1['ModifiedPeptide'] + da1['PrecursorCharge'].astype(str)
    
    # Try to apply RT normalization if available
    try:
        rt = pqp2.copy()
        ds1 = pd.merge(da1, rt, on=['ModifiedPeptide','PrecursorCharge'], how='outer')
        ds2 = ds1[ds1['ions'].isnull()==False]
        ds2 = ds2[ds2['NormalizedRetentionTime'].isnull()==False]
    except NameError:
        # If RT normalization was not performed
        ds2 = da1
    
    # FragPipe library columns
    cols = ['PrecursorMz', 'ProductMz', 'ProteinId', 
            'PeptideSequence', 'ModifiedPeptideSequence', 'PrecursorCharge', 
            'LibraryIntensity', 'NormalizedRetentionTime', 'ions']
    
    # Prepare datasets for merging
    try:
        ds3 = ds2[['PrecursorMz', 'ProductMz', 'Protein_name',
               'StrippedPeptide', 'ModifiedPeptide', 'PrecursorCharge',
               'LibraryIntensity', 'NormalizedRetentionTime', 'ions']]
        ds3.columns = cols
    except KeyError as e:
        print(f"Warning: Column mapping issue - {str(e)}")
        # Try a more flexible approach for column mapping
        required_cols = ['ModifiedPeptide', 'PrecursorCharge']
        for col in required_cols:
            if col not in ds2.columns:
                raise Exception(f"Required column '{col}' not found in SysteMHC libraries")
        
        # Map available columns and fill missing ones with NaN
        col_map = {
            'PrecursorMz': 'PrecursorMz', 
            'ProductMz': 'ProductMz',
            'Protein_name': 'ProteinId',
            'StrippedPeptide': 'PeptideSequence',
            'ModifiedPeptide': 'ModifiedPeptideSequence',
            'PrecursorCharge': 'PrecursorCharge',
            'LibraryIntensity': 'LibraryIntensity',
            'NormalizedRetentionTime': 'NormalizedRetentionTime'
        }
        
        ds3 = pd.DataFrame()
        for src_col, dst_col in col_map.items():
            if src_col in ds2.columns:
                ds3[dst_col] = ds2[src_col]
            else:
                ds3[dst_col] = np.nan
        
        ds3['ions'] = ds2['ions']
    
    sample_library3 = sample_library2[cols]
    
    # Merge libraries (exclude duplicates)
    ds4 = ds3[~ds3['ions'].isin(sample_library3['ions'])]
    lib_merge = pd.concat([sample_library3, ds4])
    
    # Save merged library
    merged_lib_path = os.path.join(output_dir, 'merged_Sample+SysteMHC_library.tsv')
    lib_merge.to_csv(merged_lib_path, sep='\t', index=False)
    print(f"Merged library saved to: {merged_lib_path}")
    
    return merged_lib_path

# Allow script to be run directly or imported as a module
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python for_fragpipe.py <sample_library_path> <output_dir> <systemhc_lib1> [systemhc_lib2 ...]")
        sys.exit(1)
    
    sample_library_path = sys.argv[1]
    output_dir = sys.argv[2]
    systemhc_lib_paths = sys.argv[3:]
    
    merge_libraries(sample_library_path, systemhc_lib_paths, output_dir)