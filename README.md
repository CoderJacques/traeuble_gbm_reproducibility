# GBM Multi-Omics Analysis Pipeline

This repository contains the code for the analysis in our [paper](https://www.biorxiv.org/content/10.64898/2026.02.26.708154v1). It contains a collection of Jupyter notebooks and scripts for the multi-omics analysis of Glioblastoma (GBM). The analysis integrates various data types, including single-cell RNA-seq, spatial transcriptomics (Visium), whole-exome sequencing (WES), and clinical data, primarily using MOFA+ (Multi-Omics Factor Analysis v2).

## Installation

To run the analyses in this repository, you need to have Conda installed. You can create the required conda environment using the provided `mofa_analysis_new.yml` file.

```bash
conda env create -f mofa_analysis_new.yml
```

After the environment is created, you can activate it with:

```bash
conda activate mofa_analysis_new
```

## File Overview

### Main Notebooks

-   **`CGGA_pp.ipynb`**: Pre-processing of data from the CGGA (Chinese Glioma Genome Atlas) dataset.
-   **`TCGA_pp.ipynb`**: Pre-processing of data from the TCGA (The Cancer Genome Atlas) dataset.
-   **`clinical.ipynb`**: Analysis of clinical data associated with the samples.
-   **`scrnaseq.ipynb`**: Analysis of single-cell RNA-sequencing data.
-   **`visium.ipynb`**: Analysis of Visium spatial transcriptomics data.
-   **`wes_analysis.ipynb`**: Analysis of Whole Exome Sequencing data.
-   **`bulk.ipynb`**: Analysis of bulk transcriptomics data.
-   **`Pathways_and_ClinicalData_Prep.ipynb`**: Preparation of pathway annotations and clinical data for integration.
-   **`cox.ipynb`**: Survival analysis using Cox proportional hazards models.
-   **`downstream.ipynb`**: Downstream analysis of the integrated multi-omics data.
-   **`ligand_receptor.ipynb`**: Ligand-receptor interaction analysis from single-cell data.
-   **`map_scpoli.ipynb`**: Cell type mapping using scPoli.
-   **`plots_mofa.ipynb`**: Generating various plots from the MOFA analysis results.
-   **`spot_deconvolution.ipynb`**: Deconvolution of spatial transcriptomics spots.
-   **`subcluster.ipynb`**: Subclustering analysis of cell populations.

### MOFA Pipeline (`mofa/`)

This directory contains the core MOFA+ pipeline scripts.

-   **`mofa/scripts/`**: Contains Jupyter notebooks and R/Python scripts that constitute the steps of the MOFA pipeline.
    -   `00_Configuration_Update.ipynb`: Updates configuration files.
    -   `00_Data_Conversion.ipynb`: Converts data into the required format for MOFA.
    -   `01_Prepare_Pseudobulk.ipynb`: Prepares pseudobulk data from single-cell data.
    -   `02_Integrate_and_Normalize_Data_Sources.ipynb`: Integrates and normalizes data from different omics sources.
    -   `03_Run_MOFA.ipynb`: Runs the main MOFA model.
    -   `04_Downstream_Factor_Analysis.ipynb`: Performs downstream analysis on the inferred MOFA factors.
    -   `06_Downstream_Pathways.ipynb`: Pathway analysis based on MOFA factors.
    -   `MS0_Libraries.r`, `MS1_Functions.py`, `MS1_Functions.r`, `MS2_Plot_Config.r`: Helper scripts with utility functions and configurations.
-   **`mofa/configurations/`**: Contains CSV files for configuring the different steps of the analysis pipeline.

### Source Code (`src/`)

-   **`src/loaders.py`**: Contains Python functions for loading and pre-processing data.

## Citation

If you use this code or data in your research, please cite our paper:

```bibtex
@article{Traeuble2026,
  title={Spatially Integrated Multi-Omics reveals the Multicellular Landscape of Progenitor-Driven Glioblastoma Progression},
  author={Traeuble, Korbinian and Traeuble, Jakob and MOSAIC consortium and Kaminski Schierle, Gabriele S and Heinig, Matthias},
  year={2026},
  doi = {10.64898/2026.02.26.708154},
  publisher={Cold Spring Harbor Laboratory}
  URL = {https://www.biorxiv.org/content/early/2026/02/28/2026.02.26.708154},
  journal={bioRxiv},
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
