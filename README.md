# DIA-Aspire
This is a computational workflow for identification and quantifiaction of DIA immunopeptidomics data.

# Overview
This workflow contains several modules:
1. Library integration: Sample library is firstly created by database search. Then the sample allele-specific library and the allele-specific libraries from SysteMHC Atlas are integrated to obtain a comprehensive one.
2. Fragment ion selection: Among the fragment ions in the integrated library, ions other than the five ion types (b-type ions (b), y-type ions (y), a-type ions (a), neutral loss ions (n), and internal ions (m)) are first removed to obtain a intermediate library, which an intensity-based filtering strategy is implemented on, obtaining top 12 abundant ions for each precursor to result in an optimized spectral library.
3. Identification and quantifiaction: DIA-NN is used to analyze the DIA immunopeptidomics data based on the optimized spectral library.

# Installation
1. Install DIA-NN. (currently, we used DIA-NN 1.8.1, but any version is supported).
2. Install Python as it is used in DIA-Aspire. For this, we recommend to create a conda environment.
   - First, run the following cmd: `conda create --name DIA-Aspire-py39 python=3.9`
   - Then, run `conda activate DIA-Aspire-py39` to activate the environment
   - Next, run `pip install numpy, pip install pandas, pip install click, pip install scikit-learn, pip install statsmodels` to install the packags needed
   - Finally, deactivate the environment by `conda deactivate`
3. Download DIA-Aspire:
   - by `git clone https://github.com/WgShaoLab/DIA-Aspire` or download the `ZIP` file.

# Usage
1. Go into the DIA-Aspire by `cd DIA-Aspire`
2. Activate the environment by `conda activate DIA-Aspire-py39` and change the ownership by `chmod 755 -R *`
3. Run DIA-Aspire by `python3 dia_aspire.py` to get the GUI of DIA-Aspire
4. Configure the DIA-NN path with the real absolute path of DIA-NN in your computer, by default it is `/usr/diann/1.8.1/diann-1.8.1`
5. Input DIA data by selecting the folder or adding files iteratively
6. Input the sample-specific library built by `Fragpipe` or `SysteMHC-pipeline`.
7. Set the absolute path of the output 
8. Selelct the HLA allele to download the allele-specific libraries from `SysteMHC Atlas`
9. Configure the parameters used by `DIA-NN`
10. Click `Run` to start the analysis. This includes retention time alignment, libraries integration, identification and quantification. And the results will be in the directory you configured before. The name of the results are all start with `lib-base-result`.

# Contact Us
For issues in using DIA-Aspire, please report to this GitHub repository.

