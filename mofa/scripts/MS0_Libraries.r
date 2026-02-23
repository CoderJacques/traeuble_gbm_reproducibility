### Set the Python environment for reticulate before loading anything
Sys.setenv(RETICULATE_PYTHON = "/home/sagemaker-user/user-default-efs/.conda/envs/mofa_analysis/bin/python")
Sys.unsetenv("PYTHONPATH")  # Avoid conflicts
Sys.setenv(PATH = paste("/home/sagemaker-user/user-default-efs/.conda/envs/mofa_analysis/bin", Sys.getenv("PATH"), sep = ":"))

### Get path where libraries have been installed
path_to_conda_env <- "/home/sagemaker-user/user-default-efs/.conda/envs/mofa_analysis"
path_to_r_libs <- paste0(path_to_conda_env, "/lib/R/library")
print(path_to_r_libs)

### Libraries for all the project scripts
library(Seurat, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(SeuratData, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
# library(SeuratDisk, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
# library(anndata, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)  # conversion hd5f to seurat file

library(data.table, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(stringr, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(dplyr, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(tidyverse, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(reshape2, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(caret, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(gridExtra, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(scales, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)

library(ggplot2, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(ggraph, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(ggpubr, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(corrplot, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(ggokabeito, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)

library(MOFA2, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
library(MOFAdata, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)

# Load reticulate and verify Python config
library(reticulate, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
print(reticulate::py_config())  # Verify the correct Python is loaded

library(DESeq2, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)
# library(muscat, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)

library(ggokabeito, quietly = TRUE, verbose = FALSE, lib.loc = path_to_r_libs)  # Loaded twice in original, keeping last instance