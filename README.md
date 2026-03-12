<p align="center" style="margin-bottom: 0px !important;">
  <img src="https://github.com/WgShaoLab/DIA-Aspire/blob/main/dia_aspire/constants/DIA-Aspire_logo.jpg" width="100" height="100">
</p>
<h1 align="center" style="margin-top: -0px; font-size: 20px">DIA-Aspire</h1>

This is a computational workflow with deep learning rescoring for identification and quantifiaction of DIA-based immunopeptidome.

# Overview
This workflow contains several modules:
1. **Library integration**: The sample allele-specific library and the allele-specific libraries from [SysteMHC Atlas](https://systemhc.sjtu.edu.cn/) are integrated to obtain a comprehensive one.
   - **Note**: For sample allele-specific library generation, the DIA data is firstly converted to pseudo-DDA by [DIA-Umpire](https://github.com/cctsou/DIA-Umpire). Then the pseudo-DDA and expirimentally acquired DDA (if available) are combined to establish sample library by database search using [SysteMHC-pipeline](https://github.com/WShaoLab/SysteMHC-pipeline) or [FragPipe](https://fragpipe.nesvilab.org/). If FragPipe is used, [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) (for **HLA-I**) or [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) (for **HLA-II**) are needed to be used to predict binding affinity. Next, the sample library is filtered by predicted binders to generate sample allele-specific library.
2. **Fragment ion selection**: Among the fragment ions in the integrated library, ions **other than** the five ion types (b-type ions (**b**), y-type ions (**y**), a-type ions (**a**), neutral loss ions (**n**), and internal ions (**m**)) are first removed to obtain a intermediate library, which an intensity-based filtering strategy is implemented on, obtaining top **12** abundant ions for each precursor to result in an optimized spectral library.
3. **Identification and quantifiaction**: [DIA-NN](https://github.com/vdemichev/DiaNN) is used to analyze the DIA immunopeptidomics data based on the optimized spectral library.
4. **Rescore**: Deep learning-based rescore is performed on the main out of DIANN.

# Installation
1. Install [DIA-NN](https://github.com/vdemichev/DiaNN). (currently, we used DIA-NN 2.0, but any version is supported).
2. Download DIA-Aspire:
   - by `git clone https://github.com/WgShaoLab/DIA-Aspire` or download the `ZIP` file.
3. Move into DIA-Aspire by `cd /path/to/DIA-Aspire`
4. Create a **[conda](https://www.anaconda.com/) environment** to install **Python** packages.
   - First, run `conda create --name dia_aspire python=3.10 -y`
   - Then, run `conda activate dia_aspire` to activate the environment (`conda deactivate` to deactivate the environment).
   - Next, run `pip install -e .` to install DIA-Aspire.
   - Finally, run `dia-aspire-gui` to open the GUI of DIA-Aspire.
   


# Usage
1. Open the GUI by `dia-aspire-gui`.
2. Configure the DIA-NN **path** with the real absolute path of DIA-NN in your computer, by default it is `/data/DIA-NN-2.0-Academia-Linux/diann-2.0/diann-linux`.
3. Input DIA data by selecting the folder or adding files iteratively.
4. Set the absolute path of the output.
5. Input the sample-specific library built by **FragPipe** or **SysteMHC-pipeline** or **MaxQuant**.
6. Selelct the HLA allele to download the allele-specific libraries from **SysteMHC Atlas**.
7. Configure the parameters used by **DIA-NN**.
8. Click to enable rescore and choose the rescore model (**DNN** or **SVM**).
9. Click `Run` to start the analysis. This includes retention time alignment, libraries integration, identification and quantification, rescore. And the results will be in the directory you configured before. The name of the results are all start with `lib-base-result`.

# How to cite
Huang, X., Gan, Z., Cui, H., Lan, T., Liu, Y., Caron, E., & Shao, W. (2023). The SysteMHC Atlas v2.0, an updated resource for mass spectrometry-based immunopeptidomics. Nucleic acids research.(https://doi.org/10.1093/nar/gkad1068)

# Contact Us
For issues in using DIA-Aspire, please report to this GitHub repository.

