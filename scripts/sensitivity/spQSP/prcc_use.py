# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 00:37:40 2020

Analyzing PRCC results

@author: Mehdi
"""
from pathlib import Path
import sys
import importlib
import matplotlib.pyplot as plt

# import PRCC function
import PRCC as prcc
import QSP_analysis as qa

#%%
import numpy as np
from pyDOE2 import lhs
# Reading data
num_samples = 20
num_param = 15
num_readout = 9

header, data = qa.read_csv('3.csv')
param_names = header[1:num_param+1]
read_names = header[num_param+1:]
lhd = data[:,1:num_param+1].astype(float)
readout = data[:,num_param+1:].astype(float)


Rho, Pval, Sig, Pval_correct = prcc.partial_corr(lhd, readout, 1e-14,Type = 'Spearman', MTC='Bonferroni')

sig_txt = np.zeros((num_param, num_readout), dtype='U8')
sig_txt[Pval_correct<5e-2] = '*'
sig_txt[Pval_correct<1e-6] = '**'
sig_txt[Pval_correct<1e-9] = '***'
param_group = ["beige"]*10 + ["khaki"]*(num_param-10)
readout_group = ["#B4B4FF"]*2+ ["mediumslateblue"]*(num_readout-2)
importlib.reload(qa)
cm = qa.cluster_map(np.transpose(Pval), read_names,param_names, 
                  (10,6), cmap="bwr",
                  annot=np.transpose(sig_txt),
                  row_colors = readout_group,
                  col_colors = param_group,
                  col_cluster=False, row_cluster=False,col_color_labels=["ABM", "QSP"],row_color_labels=["endpoint", "pretreatment"],
                  show_dendrogram = [False, False])
cm.savefig('heat.png',format='png', dpi=600,bbox_inches='tight')




