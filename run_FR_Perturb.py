#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np # for mathematical operations on large, multi-dimensional arrays and matrices
import pandas as pd # for manipulation and analysis of DataFrames and Series
import scipy # like numpy but for optimization, integration, interpolation, linear algebra, statistics, etc.
import spams # for sparse modeling, coding and dictionary learning algorithms
import scanpy # for preprocessing, visualization, clustering, trajectory inference, and differential expression analysis of (scRNA-seq) data
import time, sys, traceback, argparse #  for handling time-related operations, system-specific parameters and functions, exception handling and traceback, and command-line argument parsing
import os # for performing various operating system-related tasks, such as file and directory manipulation
import tqdm # offers a  progress bar for iterables, loops, and multiprocessing
import statsmodels.api as sma # provides classes and functions for the estimation of different statistical models and conducting statistical tests, and statistical data exploration
import statsmodels.stats as sms #  contains various statistical functions and tests for hypothesis testing, statistical analysis, and statistical modeling
import functools # for creating and manipulating closures, partial function application, etc.
from tqdm.contrib.concurrent import thread_map # for parallel mapping of a function over an iterable using threads, providing a progress bar for monitoring the progress of the execution 

# In[34]:


class Logger(object): # a class for lightweight logging
    '''
    Lightweight logging.
    TODO: replace with logging module
    '''

    def __init__(self, fh): # initialize the logger with a file handle 'fh'
        self.log_fh = open(fh, 'w') 

    def log(self, msg): # logs a message 'msg' to 'fh'
        '''
        Print to log file and stdout with a single command.
        '''
        print(msg, file=self.log_fh)
        print(msg)
        
def regress_covariates(dat, cov_mat): # regress out covariates from an expression matrix (represented by an AnnData object)
    '''
    Regress out covariates from expression matrix in place. 
    Two orders of magnitude faster than scanpy.pp.regress_out function.
    
    Parameters
    ----------
    dat: AnnData object 
    cov_mat: Cell x covariate Dataframe
    '''
    
    # Convert categorical factors to dummy variables, transforming 
    # categorical variables into a format suitable for regression analysis
    cov_mat = pd.get_dummies(cov_mat, drop_first=True) 
    
    # centering covariates, ensuring that the  intercept term in 
    # the regression model represents the overall mean of the response 
    # variable when all covariates are set to their mean values.
    cov_means = cov_mat.values.mean(axis = 0) 
    cov_mat = cov_mat.values - cov_means[np.newaxis, :] # Center covariates
    
    cov_mat = np.c_[np.ones((cov_mat.shape[0], 1)), cov_mat] # Append intercept (a column of ones to the covariate matrix)
    
    if scipy.sparse.issparse(dat.X): # checks data sparsity
        dat.X = dat.X.todense() # if sparse, make it a dense matrix to be compatible for operations, such as scipy.linalg.lstsq

    # perform linear regression to estimate the coefficients of the covariates, using the least squares method 
    lmfit = scipy.linalg.lstsq(cov_mat, dat.X, lapack_driver='gelsy')[0]
    
    # compute residuals by subtracting the predicted values (product of covariates and the coefficients 
    # from the original expression values'dat.X' to elimate the effects of covariates.
    resids = dat.X - cov_mat.dot(lmfit) 
    
    dat.X = resids # update the values after removing the covariates effects
    
def fit_skew_norm(t_star, t_nulls, side='both'):
    '''
    Compute p-values by fitting skew normal to null distribution. 
    
    Parameters
    ----------
    t_star: Test statistic
    t_null: Null statistics
    side: Which side to compute pvalues (left, right, or both)
    
    Returns
    ---------
    P-value
    '''
    
    if t_star == 0: # no difference between t-statistic and null distribution
        p = 1 # there's no significance
    else:
        fit = scipy.stats.skewnorm.fit(t_nulls) # returns the parameters of the skew normal distribution that best fit the null statistics.
        if side == 'left': # then, compute the cumulative distribution function (CDF)
            p = scipy.stats.skewnorm.cdf(t_star, *fit) # the probability of observing a value <= t_star
        elif side == 'right': # then, compute 1 - CDF at t_star
            p = 1 - scipy.stats.skewnorm.cdf(t_star, *fit) # the probability of observing a value > t_star
        elif side == 'both': # then, compute the two-tailed p-value
            p = scipy.stats.skewnorm.cdf(t_star, *fit)
            p = 2 * np.minimum(p, 1 - p) # by doubling the minimim of left & right-tailed p-values
        else:
            raise ValueError('Wrong side')
    return p

def scale_effs(B, logmeanexp, downsample_num = 25000, log_exp_baseline = 2):
    '''
    Scale effect sizes to mean expression using LOWESS (locally weighted scatterplot smoothing). 
    
    Parameters
    ----------
    B: Perturbation x gene unscaled effect size matrix
    logmeanexp: Vector of log mean expression values to scale the effects in 'B' to
    downsample_num: Number of effects used to fit curve (default is 25,000)
    log_exp_baseline: Mean effect magnitude from this log expression is taken as the value to scale to -- a reference point
    
    Returns
    ---------
    B: Perturbation x gene scaled effect size matrix
    scale_factors: Per-gene scale factors 
    '''
    
    # data fraction calculation is used for downsampling data 
    # = minumum of 1 and the ratio of num of effects to the product of the dimensions of 'B'
    data_frac = min(1, downsample_num / np.prod(B.shape))
    
    if B.shape[1] != len(logmeanexp): # ensures that the number of genes in 'B' matches the length of logmeanexp
        raise ValueError('Number of genes differs') # if not, raise a ValueError
    rand_idx = np.c_[np.random.randint(0, B.shape[0], downsample_num), # selecting a random subset of effects from 'B' for the LOWESS curve
                     np.random.randint(0, B.shape[1], downsample_num)] # np.c_ concatenates row & column indices to form pair signifying element coordinates in a single array
    to_plot = np.c_[logmeanexp[rand_idx[:,1]], np.abs(B[rand_idx[:,0],rand_idx[:,1]])] # a subset of 'logmeanexp' and corresponding absolute values of 'B' effects are selected
    to_plot = to_plot[np.where(to_plot[:,1] != 0)[0],:] # 0-effect entries are removed from the data
    to_plot[:,1] = np.log(to_plot[:,1]) # the effect sizes are transformed to the log scale  to stabilize variance and potentially normalize distributions
    fit = sma.nonparametric.lowess(to_plot[:,1], to_plot[:,0], return_sorted=False, xvals = logmeanexp) # fit a LOWESS model to the data in to_plot
    baseline = fit[min(i for i,x in enumerate(logmeanexp) if x > log_exp_baseline)] # get a baseline value from the fitted LOWESS curve at the point where logmeanexp > log_exp_baseline value
    scale_factors = np.exp(fit - baseline) # rescale the effect sizes by exponentiating the difference between the LOWESS fit and the baseline to get new scale factors
    B = B / scale_factors # scale the original effect size matrix B by dividing each gene's effect sizes by the corresponding scale factor
    return B, scale_factors # returns the scaled effect size matrix and the per-gene scale factors

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    # divide t by 60 (for seconds to minutes, then minutes to hours) and 24 (for hours to days), accumulating the remainders and quotients into a list [d, h, m, s].
    [d, h, m, s, n] = functools.reduce(lambda ll, b: divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = '' # append the following days, hours, & minutes if > 0
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

def signif(X, n):
    '''Round elements of a pandas DF X to n significant figures'''
    def func(x):
        if x == 0:
            return 0
        else:
            return round(x, n - 1 - int(np.floor(np.log10(abs(x))))) # return the rounded non-zero value
    return X.applymap(func) # apply func to each element of the DataFrame X using the applymap method




# In[6]:


os.environ['KMP_WARNINGS'] = 'off' # disable KMP (Kernel Mode-Setting Protocol) warnings to prevent irrelevant warnings from cluttering the script's output, improving user experience
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* Factorize-Recover for Perturb-seq analysis (FR-Perturb)\n" # a header to introduce the tool's name and the author
MASTHEAD += "* by Douglas Yao 2022 \n"
MASTHEAD += "*********************************************************************\n"

parser = argparse.ArgumentParser()

# define the various command-line options that the script accepts

### Required flags
parser.add_argument('--input-h5ad', default=None, type=str,
                    help='h5ad file (from the AnnData package) containing raw gene expression counts for all cells')
parser.add_argument('--input-perturbation-matrix', default=None, type=str,
                    help='Whitespace-delimited file containing a table with columns corresponding to cells and rows corresponding to perturbations. Cells containing a given perturbation should be indicated with a "1", otherwise "0".')
parser.add_argument('--control-perturbation-name', default=None, type=str,
                    help='Comma-separated list of perturbation names that represent control perturbations')
parser.add_argument('--out', default=None, type=str,
                    help="Output prefix (including directory) for effect sizes")

### Optional
parser.add_argument('--compute-pval', default=False, action='store_true',
                    help='Whether or not to compute p-values for all effect size estimates by permutation testing')           
parser.add_argument('--rank', default=20, type=int,
                    help='Hyperparameter determining the rank of the matrix during the factorize step')
parser.add_argument('--lambda1', default=0.1, type=float,
                    help='Hyperparameter determining the sparsity of the factor matrix during the factorize step of the method. Higher value = more sparse.')
parser.add_argument('--lambda2', default=10, type=float,
                    help='Hyperparameter determining the sparsity of learned effects during the recover step of the method. Higher value = more sparse.')
parser.add_argument('--covariates', default=None, type=str,
                    help='Comma-separated list of covariate names to regress out of the expression matrix (names must match the column names in the meta-data of the h5ad object)')
parser.add_argument('--guide-pooled', default=False, action='store_true',
                    help='Runs the version of FR-Perturb that assumes data is generated from guide pooling')
parser.add_argument('--cell-pooled', default=False, action='store_true',
                    help='Runs the version of FR-Perturb that assumes data is generated from cell pooling')
parser.add_argument('--num-perms', default=10000, type=int,
                    help='Number of permutations when doing permutation testing')
parser.add_argument('--fit-zero-pval', default=False, action='store_true',
                    help='Compute p-values by fitting skew-normal distribution to null distribution (allows for p-values below 1/num_perms, but significantly increases compute time)')
parser.add_argument('--multithreaded', default=False, action='store_true',
                    help='Use multithreading to fit skew-normal distributions, which can substantially reduce compute time')
parser.add_argument('--output-factor-matrices', default=False, action='store_true',
                    help='Whether or not to output the latent gene expression factor matrices in addition to the full effect sizes matrix')
parser.add_argument('--input-factorized-mat', default=None, type=str,
                    help='Rather than inputting the expression count matrix, one can specify factorized count matrices instead. The factorized matrices must be')
parser.add_argument('--cross-validate', default=None, type=int,
                    help='Whether or not to check the self-consistency of the')


# In[ ]:


if __name__ == '__main__': # ensure that the code block is executed only when the script is run directly, not when imported as a module in another script

    args = parser.parse_args()
    # args = parser.parse_args(['--input-h5ad', 'test/Simulated_seurat.h5ad',
    #                       '--input-perturbation-matrix', 'test/Simulated_perturbations.txt',
    #                       '--control-perturbation-name', 'non-targeting',
    #                       '--out', 'test/out',
    #                       '--num-perms', '100',
    #                       '--compute-pval', '--multithreaded'])
    
    #  initialize a logging object to write logs to a file named based on the output prefix provided by the user plus '.log' 
    log = Logger(args.out + '.log') # there's an ERROR here!
    
    # try:
    defaults = vars(parser.parse_args(''))
    opts = vars(args) # converts the parsed command-line arguments into a dictionary for easy access
    non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]] # create a list of arguments that the user has explicitly set (i.e., non-default arguments)
    header = MASTHEAD # appenD the following 2 lines to the header MASTHEAD
    header += "Call: \n"
    header += './run_FR_perturb.py \\\n'
    options = ['--' + x.replace('_', '-') + ' ' + str(opts[x]) + ' \\' for x in non_defaults] #  construct a list of command-line options based on non-default arguments.
    header += '\n'.join(options).replace('True', '').replace('False', '') # append the options to the header, removing the representations of boolean True and False.
    header = header[0:-1] + '\n' # trim the last backslash and adds a newline to finalize the header
    log.log(header) # log the header
    log.log('Beginning analysis at {T}'.format(T=time.ctime())) # log the start time of the analysis
    start_time = time.time() # record the start time for later calculation of total runtime

    # check for the presence of required arguments and raise errors if they're missing
    if not args.input_h5ad:
        raise ValueError('Must specify --input-h5ad')
    if not args.input_perturbation_matrix:
        raise ValueError('Must specify --input-perturbation-matrix')
    if not args.control_perturbation_name:
        raise ValueError('Must specify --control-perturbation-name')
    if not args.out:
        raise ValueError('Must specify --out')

    if args.guide_pooled and args.cell_pooled:
        raise ValueError('Only one of --guide-pooled and --cell-pooled should be set')
    elif args.cell_pooled:
        overload_type = 'droplet'
    else:
        overload_type = 'guide' # use this version by default

    log.log('Loading input data...  ')
    dat = scanpy.read_h5ad(args.input_h5ad) # load the data
    p_mat_pd = pd.read_csv(args.input_perturbation_matrix, index_col = 0, delim_whitespace=True)
    if not dat.obs.index.equals(p_mat_pd.columns): # ensure the consistency of cell names between the expression matrix and the perturbation matrix
        # breakpoint()
        raise ValueError('Cell names in perturbation matrix do not match cell names in expression matrix')
    log.log('Done')

    # center rows of p_mat_pd -- adjust if there's multiple cells in one droplet
    if overload_type == 'droplet':
        guides_per_cell = np.array(p_mat_pd.sum(axis = 0))
        p_mat_pd = p_mat_pd.divide(guides_per_cell[np.newaxis, :])

    # get covariates if specified
    if args.covariates:
        cov_names = args.covariates.split(',')
        cov_mat = dat.obs[cov_names]
        
    # normalize, log-transform, and regress out the expression matrix

    # normalization and initial processing
    log.log('Regressing out covariates and centering expression matrix...  ')
    scanpy.pp.normalize_total(dat, target_sum = 10000)
    epsilon = 1e-10 # Rawan added this
    logmeanexp = np.squeeze(np.array(np.log(np.mean(dat.X, axis = 0) + epsilon)))

    scanpy.pp.log1p(dat) # apply log transformation
    
    # egression of covariates if specified
    if args.covariates:
        regress_covariates(dat, cov_mat)
    else:
        dat.X = dat.X - dat.X.mean(axis = 0)

    # center expression matrix based on control expression
    n_guides = p_mat_pd.values.sum(axis = 0)
    nt_names = args.control_perturbation_name.split(',')
    
    # Rawan added this as we are dealing with only one control, so we've one dimension, and we can't use .values
    if(len(nt_names) == 1):
        nt_names = nt_names[0]
        ctrl_idx = np.logical_and(n_guides == 1, p_mat_pd.loc[nt_names].sum(axis = 0) != 0)
    else:
        ctrl_idx = np.logical_and(n_guides == 1, p_mat_pd.loc[nt_names].sum(axis = 0).values != 0)
    
    ctrl_exp = dat.X[ctrl_idx,:].mean(axis = 0)
    dat.X = dat.X - ctrl_exp
    log.log('Done')

    # Factorize
    log.log('Factorizing expression matrix... ')
    
    # ensure that the zeros are consistently handled across `dat.X` and `logmeanexp` -- Rawan added this, then everything worked!
    non_zero_genes = np.squeeze(np.array(dat.X.sum(axis=0))) != 0
    dat = dat[:, non_zero_genes]  # filter out zero-expression genes from dat
    logmeanexp = logmeanexp[non_zero_genes]  # make sure logmeanexp vector is filtered similarly

    keep_cells = p_mat_pd.sum(axis = 0) > 0
    p_mat = np.asfortranarray(p_mat_pd.loc[:, keep_cells].T).astype(np.float32)
    W = spams.trainDL(np.asfortranarray(dat.X.T), K=args.rank, lambda1=args.lambda1, iter=50, verbose=False)
    U_tilde = spams.lasso(np.asfortranarray(dat.X.T), D=W, lambda1=args.lambda1, verbose=False)
    U_tilde = U_tilde[:, keep_cells]
    U_tilde = np.asfortranarray(U_tilde.T.todense()).astype(np.float32)
    W = W.T
    log.log('Done')

    # Recover
    log.log('Regressing left factor matrix on perturbation design matrix...  ')
    U = spams.lasso(U_tilde, D=p_mat, lambda1=args.lambda2, verbose=False)
    B = U.dot(W)
    log.log('Done')

    # Compute pvalues by permutation testing
    if args.compute_pval:
        if not args.fit_zero_pval:
            log.log('Computing p-values by permutation testing ({} total permutations)...  '.format(args.num_perms))
            pvals = np.zeros((B.shape))
            for i in tqdm.tqdm(range(args.num_perms)):
                p_mat_perm = np.asfortranarray(p_mat[np.random.permutation(p_mat.shape[0]),:])
                U_perm = spams.lasso(U_tilde, D=p_mat_perm, lambda1=args.lambda2, verbose=False)
                B_perm = U_perm.dot(W)
                temp_indices = B < B_perm
                pvals[temp_indices] = pvals[temp_indices] + 1
            pvals /= args.num_perms
            pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5] # get 2-sided pvalues
            pvals *= 2 
            pvals = (pvals * args.num_perms + 1) / (args.num_perms + 1)
            log.log('Done')
        else:
            args.num_perms = 500
            log.log('Computing p-values by permutation testing ({} total permutations)...  '.format(args.num_perms))
            B_perms = np.empty((np.product(B.shape), args.num_perms))
            
            for i in tqdm.tqdm(range(args.num_perms)):
                p_mat_perm = np.asfortranarray(p_mat[np.random.permutation(p_mat.shape[0]),:])
                U_perm = spams.lasso(U_tilde, D=p_mat_perm, lambda1=args.lambda2, verbose=False)
                B_perm = U_perm.dot(W)
                B_perms[:,i] = np.ravel(B_perm)
            pvals = (B_perms < np.ravel(B)[:,np.newaxis]).sum(axis=1) / B_perms.shape[1]
            pvals[pvals > 0.5] = 1 - pvals[pvals > 0.5] # get 2-sided pvalues
            pvals *= 2 

            log.log('Fitting skew-normal distribution to effects with p=0 ({} total effects)...  '.format(np.sum(pvals == 0)))
            zero_indices = np.where(pvals == 0)[0]
            B_flattened = np.ravel(B)
            
            if args.multithreaded:
                fitted_pvals = thread_map(lambda i: fit_skew_norm(B_flattened[i], B_perms[i,:]), zero_indices) 
                
            else: 
                fitted_pvals = []
                for i in tqdm.tqdm(zero_indices):
                    t_star = B_flattened[i]
                    t_nulls = B_perms[i,:]
                    p = fit_skew_norm(t_star, t_nulls)
                    fitted_pvals.append(p)
            pvals[zero_indices] = fitted_pvals
            pvals = np.reshape(pvals, B.shape)
            log.log('Done')

        qvals = sms.multitest.multipletests(pvals.flatten(), method = 'fdr_bh')[1]
        qvals = np.reshape(qvals, pvals.shape)
        pvals = pd.DataFrame(data = np.transpose(pvals), index = dat.var.index, columns = p_mat_pd.index)
        qvals = pd.DataFrame(data = np.transpose(qvals), index = dat.var.index, columns = p_mat_pd.index)
        pvals = signif(pvals, 3)
        qvals = signif(qvals, 3)

    log.log('Scaling effects...  ')
    B,_ = scale_effs(B, logmeanexp)
    log.log('Done')
    log.log('Outputting results...  ')
    B = pd.DataFrame(data = np.transpose(B), index = dat.var.index, columns = p_mat_pd.index)
    # B = signif(B, 3)
    B.to_csv(args.out + '_LFCs.txt', sep = '\t')
    if args.compute_pval:
        pvals.to_csv(args.out + '_pvals.txt', sep = '\t')
        qvals.to_csv(args.out + '_qvals.txt', sep = '\t')
    log.log('All done!')
    '''       
    except Exception: # catching any exceptions
        ex_type, ex, tb = sys.exc_info()
        log.log(traceback.format_exc(ex))
        raise
    '''
    # finally: # log  the completion time and the total elapsed time
    log.log('Analysis finished at {T}'.format(T=time.ctime()))
    time_elapsed = round(time.time() - start_time, 2)
    log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))


# In[ ]:




