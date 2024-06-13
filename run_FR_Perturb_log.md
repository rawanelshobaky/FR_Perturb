# Author: Rawan Elshobaky 
#### Goal: Trying to make run_FR_Perturb.py work

# Current errors as of 3/4/24

When I run this in the command line:
    ./FR-Perturb/run_FR_Perturb.py --input-h5ad ./perturb_seq_data/GSM6858448_KO_cell_pooled.h5ad --input-perturbation-matrix ./perturb_seq_data/GSM6858448_KO_cell_pooled_perturbations.txt --compute-pval --fit-zero-pval --multithreaded --out ./cell_pooled_LFCs.txt

I get this error:
    File "./FR-Perturb/run_FR_Perturb.py", line 339, in <module>
    B_perms = np.empty((np.product(B.shape), args.num_perms))
    numpy.core._exceptions.MemoryError: Unable to allocate 32.2 GiB for an array with shape (8632200, 500) and data type float64

A MemoryError is raised as a large amount of memory (32.2 GiB) is required to create an array with the shape (8632200, 500) and data type float64. My system can't handle that, so I'll consider using smaller subsets instead of the whole dataset and see if it works.

# Update 3/27

I created subsets of each of the 4 subsets and their perturbations by including the first 50 genes only. I then ran:

conda activate FR-Perturb
./FR-Perturb/run_FR_Perturb.py --input-h5ad ./subsets/cell_pooled_genes50.h5ad --input-perturbation-matrix ./subsets/cell_pooled_perts.txt --control-perturbation-name non-targeting --out ./subsets/cell_pooled

Debugging it..

Work on debugging where the nans come from/ how to fix it

# Update 4/1

Kept adding assert statements to see where things go wrong (all values become nan)
The problem is that ctrl_idx returns all false values, so it can't compare it to the control.

After running: ctrl_idx = np.logical_and(n_guides == 1, p_mat_pd.loc[nt_names].sum(axis = 0) != 0)

I get in subsets.ipynb:
[ True False  True False False  True  True False  True  True False  True
  True  True  True False False  True False False  True False  True False
  True False False  True False  True False False False False False  True
 False  True False False False False False False  True False  True False
 False False]

I get in run_FR_Perturb.py:

AAACCCAAGACCGCCT-1    False
AAACCCAAGAGCGACT-1    False
AAACCCAAGCATTGAA-1    False
AAACCCAAGCTAATCC-1    False
AAACCCAAGCTCTGTA-1    False
AAACCCAAGTCGGCCT-1    False
AAACCCAAGTGATGGC-1    False
AAACCCAAGTGCGTCC-1    False
AAACCCAAGTTGCGAG-1    False
AAACCCACAAGAATAC-1    False
AAACCCACAAGACTGG-1    False
AAACCCACACATTCGA-1    False
AAACCCACACCGTCGA-1    False
AAACCCACACGCACCA-1    False
AAACCCACACTTGGGC-1    False
AAACCCACAGACGGAT-1    False
AAACCCACAGAGATGC-1    False
AAACCCACAGGTCTCG-1    False
AAACCCACATTGTCGA-1    False
AAACCCAGTAGGACCA-1    False
AAACCCAGTCGCTCGA-1    False
AAACCCAGTGATGAAT-1    False
AAACCCAGTGTGAATA-1    False
AAACCCAGTTTAAGGA-1    False
AAACCCAGTTTGAAAG-1    False
AAACCCATCAACTACG-1    False
AAACCCATCAGCAGAG-1    False
AAACCCATCATTGCTT-1    False
AAACCCATCCATAGAC-1    False
AAACCCATCCCGAGGT-1    False
AAACCCATCCGATTAG-1    False
AAACCCATCCTGGCTT-1    False
AAACCCATCGGTTGTA-1    False
AAACCCATCTCAGGCG-1    False
AAACCCATCTCCCAAC-1    False
AAACCCATCTCTGACC-1    False
AAACGAAAGAAGAACG-1    False
AAACGAAAGACTTAAG-1    False
AAACGAAAGAGAGTGA-1    False
AAACGAAAGCCTCAAT-1    False
AAACGAAAGCCTGGAA-1    False
AAACGAAAGCGATTCT-1    False
AAACGAAAGCGCAATG-1    False
AAACGAAAGCTCCATA-1    False
AAACGAAAGGCTGGAT-1    False
AAACGAAAGGTAAGTT-1    False
AAACGAAAGGTTCTTG-1    False
AAACGAAAGTAGCTCT-1    False
AAACGAAAGTGCAACG-1    False
AAACGAAAGTGCGACA-1    False
dtype: bool

I confirmed that I have the same n_guides, nt_names, and p_mat_pd, so I'm not sure where the discrepancy arises from.

Fixed it! By picking the right control name (nt_names)

There's now a new error:

File "./FR-Perturb/run_FR_Perturb.py", line 321, in <module>
    dat.X = dat.X[:,np.squeeze(np.array(dat.X.sum(axis = 0))) != 0]
  File "/home/ec2-user/miniforge3/envs/FR-Perturb/lib/python3.8/site-packages/anndata/_core/anndata.py", line 692, in X
    raise ValueError(
ValueError: Data matrix has wrong shape (50, 11494), need to be (50, 14387).

# Update 4/10

The shape of dat.X doesn't change anywhere except in squeeze when the columns that only have 0's are removed.
- This throws an error. Am I supposed to not have any columns with just 0's?
- To get around that, I tried changing the 0's to nan (I'm not sure if that's a good idea), but that didn't work either:(

Fixed it by updating dat instead of dat.X
Now, have to fix the issue of different magnitudes in the scaling function.

# Update 4/16

Fixed the issue by removing the zero genes from both dat (perturbation x gene unscaled effect size matrix) and logmeanexp (log mean expression values vector) so that they could be scaled together later after factorization

It now fully works!!

The output [INPUT]_LFCs.txt is a whitespace-delimited text file that shows the log-fold changes in gene expression compared to the expression levels in cells with control guides where genes correspond to rows and perturbations correspond to columns.

I ran the code for the four subsets I made:
  [1] ./FR-Perturb/run_FR_Perturb.py --input-h5ad ./subsets/cell_pooled_genes50.h5ad --input-perturbation-matrix ./subsets/cell_pooled_perts.txt --control-perturbation-name non-targeting --out ./subsets/cell_pooled
  [2] ./FR-Perturb/run_FR_Perturb.py --input-h5ad ./subsets/conventional_cp_genes50.h5ad --input-perturbation-matrix ./subsets/conventional_cp_perts.txt --control-perturbation-name non-targeting --out ./subsets/conventional_cp
  [3] ./FR-Perturb/run_FR_Perturb.py --input-h5ad ./subsets/conventional_gp_genes50.h5ad --input-perturbation-matrix ./subsets/conventional_gp_perts.txt --control-perturbation-name non-targeting --out ./subsets/conventional_gp
  [4] ./FR-Perturb/run_FR_Perturb.py --input-h5ad ./subsets/guide_pooled_genes50.h5ad --input-perturbation-matrix ./subsets/guide_pooled_perts.txt --control-perturbation-name non-targeting --out ./subsets/guide_pooled

There were problems with the other subsets, so I have to fix my subsets for that first

Next step, plug into Eric's benchmarking software (sent on slack): https://github.com/ekernf01/ggrn/tree/main/ggrn_docker_backend
  - write a python script that can read the input in the expected format and output it in the expected format
  - use the ggrn conda environment
  - make a directory named fake_docker
  - expected data format: check what ggrn.validate_training_data(train) checks for.
    find that in ggrn/docs/reference.md

# Update 4/24

The problem was with the subsets.
- Perturbations were correctly subset to only have 50 genes.
- 50 genes subset from the h5ad wasn't subset correctly, but I fixed that.
- It works for cell_pooled only, not for anything else that end up having not correlating values

Write comments for Eric's code (sent on slack): https://github.com/ekernf01/pereggrn_perturbations/blob/69a007eaa87bde18128bf24c743c8ffed071084e/pereggrn_perturbations.py#L41C5-L41C31

Schedule a meeting for Alexis and Eric to update on what I've been working on this semester:




conda activate FR-Perturb
