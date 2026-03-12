# maxquant_api.py - API for handling MaxQuant generated libraries

import re
import operator
import numpy as np
import pandas as pd
import os
import sys
import dia_aspire.converters.msms2tsv as msms2tsv
import dia_aspire.converters.irt_alignment as irt_align

def merge_libraries(sample_library_path, systemhc_lib_paths, output_dir):
    """
    Merge sample library with SysteMHC libraries for MaxQuant pipeline
    
    Parameters:
    -----------
    sample_library_path : str
        Path to the sample library MSMS file (msms.txt)
    systemhc_lib_paths : list
        List of paths to SysteMHC library TSV files
    output_dir : str
        Directory where the merged library will be saved
    
    Returns:
    --------
    str
        Path to the merged library file
    """
    # print(f"Sample library: {sample_library_path}")
    # print(f"SysteMHC libraries: {systemhc_lib_paths}")
    # print(f"Output directory: {output_dir}")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load sample library (MaxQuant msms.txt format)
    try:
        sample_library = msms2tsv.convert_msms2tsv(sample_library_path)
        sample_library2 = sample_library.copy()
        sample_library2['ions'] = sample_library2['ModifiedPeptide'] + sample_library2['PrecursorCharge'].astype(str)
        sample_library2['NormalizedRetentionTime'] = sample_library2['iRT']
        sample_convert_path = os.path.join(output_dir, 'Sample_library_msmstxt2tsv.tsv')
        sample_library2.to_csv(sample_convert_path,sep='\t',index=False)
        # print(f"Successfully loaded sample library: {sample_library_path}")
    except Exception as e:
        raise Exception(f"Failed to load sample library: {str(e)}")
    
    # Load irt data for RT normalization (if available)
    try:
        # Try to use the irt file in the output directory
        irt_file_path = os.path.join(output_dir, 'irt_SYSTEMHC.parquet')
        if os.path.exists(irt_file_path):
            df_need = pd.read_parquet(irt_file_path)
        else:
            # if file doesn't exist in output dir
            print(f"Failed:  {irt_file_path} not found \n")
            raise FileNotFoundError(f"{irt_file_path} not found \n")
        
        # RT normalization
        rt_reference_run = sample_library2[['ModifiedPeptide','PrecursorCharge','NormalizedRetentionTime']].drop_duplicates()
        rt_reference_run.columns = ['modified_peptide','precursor_charge','irt']
    
        df_need.columns = ['modified_peptide','precursor_charge','RT']
    
        aligned_runs1 = irt_align.lowess2(df_need, rt_reference_run, 'RT', 'irt', 0.01, 0, 10)
        pepida1 = aligned_runs1
        pepida1 = pepida1.loc[np.isfinite(pepida1['irt'])]
        pqp = pepida1
        pqp2 = pqp.groupby(['modified_peptide','precursor_charge'])[['irt']].median().reset_index()
        pqp2.columns = ['ModifiedPeptide','PrecursorCharge','NormalizedRetentionTime']
    
        # Save RT alignment results
        rt_aligned_path = os.path.join(output_dir, 'rt_aligned2reference_maxquant.csv')
        pqp2.to_csv(rt_aligned_path, index=False)
        # print(f"RT alignment results saved to: {rt_aligned_path}")
    except Exception as e:
        print(f"Warning: RT normalization skipped - {str(e)} \n")
        # Continue without RT normalization
        pqp2 = None
    
    # Combine SysteMHC libraries
    systemhc_libs = []
    for libp in systemhc_lib_paths:
        try:
            datmp = pd.read_csv(libp, sep='\t')
            systemhc_libs.append(datmp)
            # print(f"Loaded SysteMHC library: {libp}")
        except Exception as e:
            print(f"Warning: Failed to load {libp} - {str(e)} \n")
    
    if not systemhc_libs:
        raise Exception("No valid SysteMHC libraries were loaded \n")
    
    da = pd.concat(systemhc_libs)
    da = da.drop_duplicates()
    da1 = da.copy()
    da1['ions'] = da1['ModifiedPeptide'] + da1['PrecursorCharge'].astype(str)
    
    # Try to apply RT normalization if available
    try:
        if pqp2 is not None:
            rt = pqp2.copy()
            ds1 = pd.merge(da1, rt, on=['ModifiedPeptide','PrecursorCharge'], how='outer')
            ds2 = ds1[ds1['ions'].isnull()==False]
            ds2 = ds2[ds2['NormalizedRetentionTime'].isnull()==False]
        else:
            ds2 = da1
    except Exception as e:
        print(f"Warning: Error during RT normalization application - {str(e)} \n")
        ds2 = da1
    
    # MaxQuant pipeline columns
    cols = ['PrecursorMz', 'ProductMz', 'uniprot_id', \
            'StrippedPeptide', 'ModifiedPeptide', 'PrecursorCharge', 'FragmentCharge','FragmentType','FragmentNumber',\
            'LibraryIntensity', 'NormalizedRetentionTime', 'shared', 'decoy', 'ions']
    
    # Prepare datasets for merging
    try:
        ds3 = ds2[cols]
    except KeyError as e:
        print(f"Warning: Column missing in SysteMHC libraries - {str(e)} \n")
        # For any missing columns, add them with default values
        for col in cols:
            if col not in ds2.columns:
                if col == 'shared':
                    ds2[col] = 0
                elif col == 'decoy':
                    ds2[col] = 0
                else:
                    ds2[col] = np.nan
        ds3 = ds2[cols]
    
    sample_library3 = sample_library2[cols]
    
    # Merge libraries (exclude duplicates)
    ds4 = ds3[~ds3['ions'].isin(sample_library3['ions'])]
    lib_merge = pd.concat([sample_library3, ds4])
    
    # Save merged library
    merged_lib_path = os.path.join(output_dir, 'merged_Sample+SysteMHC_library_maxquant.tsv')
    lib_merge.to_csv(merged_lib_path, sep='\t', index=False)
    # print(f"Merged library saved to: {merged_lib_path}")
    
    return merged_lib_path

# Allow script to be run directly or imported as a module
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python maxquant_api.py <sample_library_path> <output_dir> <systemhc_lib1> [systemhc_lib2 ...] \n")
        sys.exit(1)
    
    sample_library_path = sys.argv[1]
    output_dir = sys.argv[2]
    systemhc_lib_paths = sys.argv[3:]
    
    merge_libraries(sample_library_path, systemhc_lib_paths, output_dir)
