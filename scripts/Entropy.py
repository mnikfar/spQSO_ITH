import os
import pandas as pd
import numpy as np
import glob2 as glob
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
result_dir = './snapShots/'
patient_spatial = pd.read_csv(result_dir+'score_840.csv')
patient_slice = patient_spatial.loc[patient_spatial['y'] == 5]
immune_cell = patient_slice.loc[patient_slice['Type'] != 1]
cancer_cell = patient_slice.loc[patient_slice['Type'] == 1]
p = len(cancer_cell)/(len(cancer_cell)+len(immune_cell))
xc=[]
zc=[]
for i in range(len(cancer_cell)):
    xc.append(cancer_cell.iloc[i]['x']) 
    zc.append(cancer_cell.iloc[i]['z'])
xi=[]
zi=[]    
for i in range(len(immune_cell)):
    xi.append(immune_cell.iloc[i]['x']) 
    zi.append(immune_cell.iloc[i]['z']) 
E = 0
sin = 0
sex = 0
for i in range(len(cancer_cell)):
    xin=(xc-xc[i])**2
    yin = (zc-zc[i])**2
    din=(xin+yin)**(0.5)
    sin=sin+din.sum()*1e-6
    xex=(xi-xc[i])**2
    yex = (zi-zc[i])**2
    dex=(xex+yex)**(0.5)
    sex=sex+dex.sum()*1e-6
E=E-(sin/sex)*p*np.log(p)
print(E)