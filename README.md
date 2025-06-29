<p align="center" style="margin-bottom: 0px !important;">
  <img src="https://github.com/WgShaoLab/DIA-Aspire/blob/master/src/DIA-Aspire_logo.jpg" width="100" height="100">
</p>
<h1 align="center" style="margin-top: -0px; font-size: 20px">DIA-Aspire</h1>

This is a computational workflow for identification and quantifiaction of DIA-based immunopeptidome.

# Overview
This workflow contains several modules:
1. Library integration: The sample allele-specific library and the allele-specific libraries from [SysteMHC Atlas](https://systemhc.sjtu.edu.cn/) are integrated to obtain a comprehensive one.
   - **Note**: For sample allele-specific library generation, the DIA data is firstly converted to pseudo-DDA by [DIA-Umpire](https://github.com/cctsou/DIA-Umpire). Then the pseudo-DDA and expirimentally acquired DDA (if available) are combined to establish sample library by database search using [SysteMHC-pipeline](https://github.com/WShaoLab/SysteMHC-pipeline) or [FragPipe](https://fragpipe.nesvilab.org/). If FragPipe is used, [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) (for **HAL-I**) or [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) (for **HLA-II**) are needed to be used to predict binding affinity. Next, the sample library is filtered by predicted binders to generate sample allele-specific library.
3. Fragment ion selection: Among the fragment ions in the integrated library, ions **other than** the five ion types (b-type ions (**b**), y-type ions (**y**), a-type ions (**a**), neutral loss ions (**n**), and internal ions (**m**)) are first removed to obtain a intermediate library, which an intensity-based filtering strategy is implemented on, obtaining top **12** abundant ions for each precursor to result in an optimized spectral library.
4. Identification and quantifiaction: [DIA-NN](https://github.com/vdemichev/DiaNN) is used to analyze the DIA immunopeptidomics data based on the optimized spectral library.

# Installation
1. Install [DIA-NN](https://github.com/vdemichev/DiaNN). (currently, we used DIA-NN 1.8.1, but any version is supported).
2. Install **Python** as it is used in DIA-Aspire. For this, we recommend to create a **[conda](https://www.anaconda.com/) environment**.
   - First, run the following cmd: `conda create --name DIA-Aspire-py39 python=3.9`
   - Then, run `conda activate DIA-Aspire-py39` to activate the environment
   - Next, run `pip install numpy, pip install pandas, pip install click, pip install scikit-learn, pip install statsmodels` to install the packags needed
   - Finally, deactivate the environment by `conda deactivate`
3. Download DIA-Aspire:
   - by `git clone https://github.com/WgShaoLab/DIA-Aspire` or download the `ZIP` file.

# Usage
1. Unzip the file of SysteMHC retention time by `cd DIA-Aspire/src` and `unzip irt_SYSTEMHC.zip`
2. Then, move into the DIA-Aspire by `cd ../`
3. Activate the environment by `conda activate DIA-Aspire-py39` and change the ownership by `chmod 755 -R *`
4. Run DIA-Aspire by `python3 dia_aspire.py` to get the GUI of DIA-Aspire
5. Configure the DIA-NN **path** with the real absolute path of DIA-NN in your computer, by default it is `/usr/diann/1.8.1/diann-1.8.1`
6. Input DIA data by selecting the folder or adding files iteratively
7. Input the sample-specific library built by **FragPipe** or **SysteMHC-pipeline**.
8. Set the absolute path of the output 
9. Selelct the HLA allele to download the allele-specific libraries from **SysteMHC Atlas**
10. Configure the parameters used by **DIA-NN**
11. Click `Run` to start the analysis. This includes retention time alignment, libraries integration, identification and quantification. And the results will be in the directory you configured before. The name of the results are all start with `lib-base-result`.

# How to cite
Huang, X., Gan, Z., Cui, H., Lan, T., Liu, Y., Caron, E., & Shao, W. (2023). The SysteMHC Atlas v2.0, an updated resource for mass spectrometry-based immunopeptidomics. Nucleic acids research.(https://doi.org/10.1093/nar/gkad1068)

# Contact Us
For issues in using DIA-Aspire, please report to this GitHub repository.

