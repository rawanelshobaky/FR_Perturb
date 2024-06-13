# Author: Rawan Elshobaky 
- I summarized the work from https://www.nature.com/articles/s41587-023-01964-9. I debugged and documented their run_FR_Perturb.py code.

# Date: 3/11/24

# Factorize-Recover for Perturb-seq analysis (FR-Perturb)

FR-Perturb is a tool operated via the command line designed to estimate the effect sizes on gene expression from Perturb-seq datasets, utilizing the factorize-recover approach [by Sharan and et al. 2019 ICML] (http://proceedings.mlr.press/v97/sharan19a/sharan19a-supp.pdf)

Perturb-seq is useful for the systematic introduction of genetic perturbations to cells and detailed analysis of the resulting changes in gene expression at the single-cell level. However, it has high costs, limitations in cell availability, and the impracticality of exhaustively screening for genetic interactions due to the vast number of genes and potential interactions. Instead, the inference method FR-Perturb can be used to  learn individual perturbation effects from composite samples through compressed Perturb-seq as outlined below.

This preliminary code is ready to use, but there's an error as it runs out of memory when Logger is called -- because some machines do not have the sufficient memory to take all the data, so you can subset them.

# Understanding the Data
Three types of data are considered:
- They all have the same `vars`: '_index', 'features'
- They all have the same `obs` except for the second conventional dataset and the guide-pooled as they have an extra `obs` ('Biological_replicate'). For example:
    - "GSM6858448_KO_cell_pooled.h5ad" is an AnnData object with n_obs × n_vars = 86956 × 16952
    obs: 'Total_RNA_count', 'Total_unique_genes', '10X_channel', 'Percent_mitochondrial_reads', 'Guides', 'Guides_collapsed_by_gene', 'Total_number_of_guides', 'S_score', 'G2M_score', 'Cell_cycle_phase'
    - "GSM6858450_KD_guide_pooled.h5ad" is an AnnData object with n_obs × n_vars = 24192 × 15668
    obs: 'Biological_replicate', 'Total_RNA_count', 'Total_unique_genes', 'Percent_mitochondrial_reads', 'Guides', 'Guides_collapsed_by_gene', 'Total_number_of_guides', '10X_channel', 'S_score', 'G2M_score', 'Cell_cycle_phase'

1. Conventional Perturb-seq: 
    - It's a traditional approach where each cell in a pooled CRISPR screen is perturbed with a single guide RNA (sgRNA) that targets a specific gene. The changes in gene expression are then analyzed using single-cell RNA sequencing (scRNA-seq)
    - This is used as a baseline for comparing the efficiency of the other new approaches introduced below.
    - There are two used datasets that have different dimensions and `obs` to be used for comparing cell_pooled and guide_pooled individually:
        - "GSM6858447_KO_conventional.h5ad" is an AnnData object with n_obs × n_vars = 86956 × 16952 
        obs: 'Total_RNA_count', 'Total_unique_genes', '10X_channel', 'Percent_mitochondrial_reads', 'Guides', 'Guides_collapsed_by_gene', 'Total_number_of_guides', 'S_score', 'G2M_score', 'Cell_cycle_phase'
        - "GSM6858449_KD_conventional.h5ad" is an AnnData object with n_obs × n_vars = 66283 × 18017
        obs: 'Biological_replicate', 'Total_RNA_count', 'Total_unique_genes', 'Percent_mitochondrial_reads', 'Guides', 'Guides_collapsed_by_gene', 'Total_number_of_guides', '10X_channel', 'S_score', 'G2M_score', 'Cell_cycle_phase'
2. Cell-pooled Perturb-seq:
    - Multiple cells are pooled within the same droplet during the sequencing process to capture a composite signal from multiple perturbations in a single readout. The signal is then computationally decompressed to estimate the effects of individual perturbations.
    - This method utilizes the sparsity (inherent simplicity where a limited number of elements control the majority of genes) of cellular regulatory networks . Consequently, it improves the efficiency while reducing the cost of the screenings.
3. Guide-pooled Perturb-seq:
    - Multiple sgRNAs are pooled within a single cell, leading to potentially introducing multiple perturbations per cell. Computational methods are then used to separate and identify the individual effects.
    - It also improves the efficiency and cost by utilizing the sparsity of cellular regulatory networks but at the level of perturbations within each cell rather than pooling cells.

# Code Explanation

The following classes/functions are defined:
- regress_covariates to regress covariates through linear regression and the least squares method 
- fit_skew_norm to compute p-values by fitting skew normal to null distribution using test statistic and null statistics
- scale_effs to scale effect sizes to mean expression using LOWESS (locally weighted scatterplot smoothing) by comparing a subset of 'B' (perturbation x gene unscaled effect size matrix) and a vector of log mean expression values to scale the effects in 'B' to, then fitting a LOWESS model
- Logger, sec_to_str, and signif to log the process and keep track of the timing

When the program runs, it:
- identifies the required flags from the user input on the command line:
    - `--input-h5ad` should be a h5ad file (from the AnnData package) containing raw gene expression counts for all cells
    - `--input-perturbation-matrix` should be a whitespace-delimited file containing a table with columns corresponding to cells and rows corresponding to perturbations. Cells containing a given perturbation should be indicated with a "1", otherwise "0".
    - `--control-perturbation-name` should be a comma-separated list of perturbation names that represent control perturbations
    - `--out` should be an output prefix (including directory) for effect sizes 

- identifies the optional flags, if any, from the user input on the command line:
    - `--compute-pval` should be whether or not to compute p-values for all effect size estimates by permutation testing
    - `--rank` should be a hyperparameter determining the rank of the matrix during the factorize step
    - `--lambda1` should be a hyperparameter determining the sparsity of the factor matrix during the factorize step of the method. Higher value = more sparse
    - `--lambda2` should be a hyperparameter determining the sparsity of learned effects during the recover step of the method. Higher value = more sparse
    - `--covariates` should be a comma-separated list of covariate names to regress out of the expression matrix (names must match the column names in the meta-data of the h5ad object)
    - `--guide-pooled` should be run the version of FR-Perturb that assumes data is generated from guide pooling 
    - `--cell-pooling` should be run the version of FR-Perturb that assumes data is generated from cell pooling 
    - `--num-perms` should be the number of permutations when doing permutation testing 
    - `--fit-zero-pval` should be a boolean to compute p-values by fitting skew-normal distribution to null distribution (allows for p-values below 1/num_perms, but significantly increases compute time)
    - `--multithreaded` should be a boolean to use multithreading to fit skew-normal distributions, which can substantially reduce compute time

- The following optional flags were added in the code, yet the author doesn't mention them in the README or use them through the rest of the code
    - `--output-factor-matrices` should be whether or not to output the latent gene expression factor matrices in addition to the full effect sizes matrix
    - `--input-factorized-mat` should be a matrix to specify factorized count matrices instead of inputting the expression count matrix
    - `--cross-validate`

- checks that all required flagss are present, and throws a ValueError if not
- loads the data from `input_h5ad` and `input_perturbation_matrix` (the expression matrix), ensure the consistency of cell names between them, and throw an Error if they're inconsistent
- adjusts the expression matrix if there's multiple cells in one droplet by centering its rows
- regresses out covariates using the regress_covariates function
- adjusts the the expression matrix to account for baseline expression levels by centering the data around the expression profile of control cells
- factorizes the expression matrix into two matrices (W and U_tilde) through dictionary learning and LASSO regression
- estimates the specific effects of genetic perturbations by performing a regression analysis using the LASSO (Least Absolute Shrinkage and Selection Operator) method
- computes p-values by permutation testing based on the value of `fit_zero_pval`, fits skew-normal distribution, then scales the effects
- logs the completion time and the total elapsed time

# Setup, Usage, & Expected Output

To set up the required environment:
- clone the repository:
    git clone https://github.com/douglasyao/FR-Perturb.git
    cd FR-Perturb
- activate the conda environment
    conda env create --file environment.yml
    conda activate FR-Perturb

Then, run the following command on the command line to estimate effect sizes using the parameters explained above:
    ./run_FR_Perturb.py --input-h5ad [INPUT_H5AD] --input-perturbation-matrix [INPUT_PERTURBATION_MATRIX] --control-perturbation-name [CONTROL_PERTURBATION_NAME] --covariates [COVARIATES] --compute-pval --fit-zero-pval --multithreaded --out [OUT]

The expected output should be a whitespace-delimited text file that shows the log-fold changes in gene expression compared to the expression levels in cells with control guides where genes correspond to rows and perturbations correspond to columns. The output name will be the prefix with _LFCs.txt appended to it. When `--compute-pval` is true, FR-Perturb also produces a file with p-values corresponding to the same gene-perturbation layout as the effect size file, as well as a file with q-values calculated through the Benjamini-Hochberg method. Note the following example:

```
PERT_1  PERT_2  PERT_3
GENE_1  0.39 -0.9 2.1 
GENE_2  0.32 -0.82 3.1
GENE_3  -0.61 0.24 1.2
```
