# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 00:37:40 2020

Analyzing PRCC results

@author: Chang
"""
from pathlib import Path
import sys
import importlib
import matplotlib.pyplot as plt

# import PRCC function
path_analysis = Path('D:/Research/Source/Repo/PopelServer/scripts/Python/analysis')
sys.path.append(str(path_analysis))
import PRCC as prcc
import QSP_analysis as qa

#%%
import numpy as np
from pyDOE2 import lhs
# generate data

num_samples = 400
num_param = 20
num_readout = 5

param_names = ['param_{}'.format(x) for x in range(num_param)]
read_names = ['readout_{}'.format(x) for x in range(num_readout)]
lhd = lhs(num_param, samples=num_samples)
readout = np.random.rand(num_samples, num_readout)

readout[:,2] = lhd[:,4] + readout[:,2]/10
readout[:,4] = -lhd[:,10]/2 + readout[:,4]
readout[:,3] = -lhd[:,9] + lhd[:,8]/5 + readout[:,3]
#%% PRCC

Rho, Pval, Sig, Pval_correct = prcc.partial_corr(lhd, readout, 1e-3,Type = 'Spearman', MTC='Bonferroni')
#%%
sig_txt = np.zeros((num_param, num_readout), dtype='U8')
sig_txt[Pval_correct<5e-2] = '*'
sig_txt[Pval_correct<1e-6] = '**'
sig_txt[Pval_correct<1e-9] = '***'

param_group = ["beige"]*5 + ["khaki"]*(num_param-5)
readout_group = ["#B4B4FF"]*2 + ["#E6E6FF"] + ["mediumslateblue"]*(num_readout-3)

#%%
importlib.reload(qa)
cm = qa.cluster_map(Rho, param_names, read_names, 
                 (5,6), cmap="RdBu",
                 annot=sig_txt,
                 row_colors = param_group,
                 col_colors = readout_group,
                 show_dendrogram = [True, True])
cm.savefig(data_folder+'heat.svg',format='svg', dpi=DPI,bbox_inches='tight')
