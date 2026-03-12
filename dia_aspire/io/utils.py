from pathlib import Path
import pandas as pd

def find_ms_files(directory: str, file_type: str) -> dict[str, str]:
    """Find MS files in directory by extension. Returns {raw_name: file_path}."""
    ext_map = {"mzml": [".mzML", ".mzml"], "hdf5": [".hdf5", ".hdf"]}
    extensions = ext_map.get(file_type, [f".{file_type}"])
    files = []
    for ext in extensions:
        pattern = f"*{ext}"
        for file in Path(directory).glob(pattern):
            if file.is_file():
                files.append(str(file))
    return ms_files_to_dict(files)


def ms_files_to_dict(ms_files: list[str]) -> dict[str, str]:
    """Convert list of MS file paths to {raw_name: file_path} dict."""
    result = {}
    for f in ms_files:
        p = Path(f)
        # Remove all extensions (e.g., "sample.mzML.hdf5" -> "sample")
        raw_name = p.name
        while p.suffix:
            raw_name = p.stem
            p = Path(raw_name)
        result[raw_name] = f
    return result

def get_filtered(data: pd.DataFrame, fdrv: float = 0.01, MBR: bool = False):
    df = data.copy()
    if MBR:
        df = df[df['fdr2_search2']<=fdrv] # Lib.PG.Q.Value
        df = df[df['fdr3_search1']<=fdrv] # Q.Value
    else:
        df = df[df['fdr1_search1']<=fdrv]  # Global.Q.Value
        df = df[df['fdr2_search1']<=fdrv] # Global.PG.Q.Value
        df = df[df['fdr3_search1']<=fdrv] # Q.Value   
    return df

def correct_fdr(data: pd.DataFrame, fdrv: float = 0.01, MBR: bool = False):
    data1 = get_filtered(data,fdrv,MBR)
    # fdr<=0.01
    data2 = data1.copy()
    data2['fdr'] = data1[['fdr1_search1','fdr2_search1','fdr3_search1']].min(axis=1)
    # fdr>0.01
    datatmp = pd.concat([data,data1])
    data3 = datatmp.drop_duplicates(keep=False)
    data3 = data3.copy()
    data3['fdr'] = data3[['fdr1_search1','fdr2_search1','fdr3_search1']].max(axis=1)
    # concat
    df = pd.concat([data2,data3])
    
    return df
