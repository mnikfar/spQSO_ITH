# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 07:56:49 2020
@author: Chang Gong

Adapted from Matlab code published in the following paper:
Marino, Simeone, Ian B. Hogue, Christian J. Ray, and Denise E. Kirschner. 
"A methodology for performing global uncertainty and sensitivity analysis 
in systems biology." Journal of theoretical biology 254, no. 1 (2008): 178-196.

"""

import numpy as np
import matlab.engine

"""
Calculate partial correlation between parameters and results
n simulations, k parameters, m results
LHS: parameter values, n x k
res: simulation results, n x m
alpha: significance level
Type = ['Spearman' (default) |'Pearson']: rank/linear correlation 
MTC: multiple testing correction, [None (default)|'Bonferroni']

output:
Rho: partial correlation, k x m
rho[i, j]: partial correlation between Y[:, j] and LHS[:, i], controlling for LHS[:, [:i,i+1:]]
Pval: p-Values, k x m
Sig: significant, k x m bool type
Pval_corrected: after multiple testing correction
"""
def partial_corr(LHS, res, alpha, Type = 'Spearman', MTC = None):

    if Type not in ['Spearman', 'Pearson']:
        raise ValueError('mode should be either "Spearman" for rank correlation or "Pearson" for linear correlation')
    
    n, k = LHS.shape
    n2 = res.shape[0]
    if n2 != n:
        raise ValueError('X and Y should have the same number of rows')
    
    if len(res.shape) == 1:
        res = res.reshape((-1, 1))
    m = res.shape[1]
    
    Rho = np.zeros((k, m))
    Pval = np.zeros((k, m))
    
    eng = matlab.engine.start_matlab()
    
    x = matlab.double(res.tolist())

    for i in range(k):
        idx = np.asarray(k * [True])
        idx[i] = False
        y = matlab.double(LHS[:, [i]].tolist())
        z = matlab.double(LHS[:, idx].tolist())
        rho, pval = eng.partialcorr(x, y, z, 'Type', Type, 'Rows','complete', nargout = 2)
        #print(rho)
        Rho[i, :] = np.array(rho).T
        Pval[i, :] = np.array(pval).T
    
    eng.quit()
    
    # Multiple testing correction: Bonferroni and FDR
    Sig = np.zeros((k, m), dtype=bool)
    Pval_correct = Pval
    
    if MTC is None:
        Sig = Pval < alpha
    elif MTC == 'Bonferroni':
        Pval_correct = Pval * m * k
        Sig = Pval_correct  < alpha 
    else:
        raise ValueError('Unknown MTC method: ' + MTC)
    
    
    return Rho, Pval, Sig, Pval_correct
#%%

# example: using matlab hospital data
if (__name__ == '__main__'):
    
    # "hospital" dataset from matlab
    eng = matlab.engine.start_matlab()
    hospital  = eng.load('hospital')['hospital']
    sexID = np.array(eng.grp2idx(eng.getfield(hospital, 'Sex')))
    age = np.array(eng.getfield(hospital, 'Age'))
    smoker = np.array(eng.getfield(hospital, 'Smoker'))
    weight = np.array(eng.getfield(hospital, 'Weight'))
    bloodpressure = np.array(eng.getfield(hospital, 'BloodPressure'))
    eng.quit()
    
    # PRCC
    X = bloodpressure
    Y = np.concatenate((weight, age, sexID, smoker), axis = 1)
    Rho, Pval, Sig, Pval_correct = partial_corr(Y, X, 1e-14, MTC='Bonferroni')
    print(np.array(Rho))
    
    