import os
import pandas as pd
import numpy as np
import glob2 as glob
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
result_dir = './snapShots/'
patient_spatial = pd.read_csv(result_dir+'score_360.csv')
patient_mixing_score = []
for y in range(20,30,1):
    mix = 0
    patient_slice = patient_spatial.loc[patient_spatial['y'] == y]
    immune_cell = patient_slice.loc[patient_slice['Type'] != 1]
    for i in range(len(immune_cell)):
        x = immune_cell.iloc[i]['x']
        z = immune_cell.iloc[i]['z']
        cell_1 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z+1)]
        cell_2 = patient_slice.loc[(patient_slice['x'] == x+1) & (patient_slice['z'] == z)]
        cell_3 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z-1)]
        cell_4 = patient_slice.loc[(patient_slice['x'] == x-1) & (patient_slice['z'] == z)]
            
        pos1 = cell_1['Type'].isin([1]).any()
        pos2 = cell_2['Type'].isin([1]).any()
        pos3 = cell_3['Type'].isin([1]).any()
        pos4 = cell_4['Type'].isin([1]).any()
            
        if (pos1 or pos2 or pos3 or pos4):
            mix += 1
            
    mixing_score = mix/(len(immune_cell) + 1) ## immune cell next to cancer cell / total immune cell
    patient_mixing_score.append(mixing_score)
print(np.mean(patient_mixing_score))