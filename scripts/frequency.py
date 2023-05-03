import os
import pandas as pd
import numpy as np
import glob2 as glob
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
result_dir = './snapShots/'
patient_spatial = pd.read_csv(result_dir+'score_240.csv')
density = []
for y in range(5,6,1):
    patient_slice = patient_spatial.loc[patient_spatial['y'] == y]
    immune_cell = patient_slice.loc[patient_slice['Type'] != 1]
    cancer_cell = patient_slice.loc[patient_slice['Type'] == 1]
    for i in range(len(immune_cell)):
        x = immune_cell.iloc[i]['x']
        z = immune_cell.iloc[i]['z']
        cell_0 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z)]
        cell_1 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z+1)]
        cell_2 = patient_slice.loc[(patient_slice['x'] == x+1) & (patient_slice['z'] == z)]
        cell_3 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z-1)]
        cell_4 = patient_slice.loc[(patient_slice['x'] == x-1) & (patient_slice['z'] == z)]
        d_0= cell_0.loc[cell_0['Type'] == 1]
        d_1= cell_1.loc[cell_1['Type'] == 1]
        d_2= cell_2.loc[cell_2['Type'] == 1]
        d_3= cell_3.loc[cell_3['Type'] == 1]
        d_4= cell_4.loc[cell_4['Type'] == 1]
        mixing_score = len(d_1)+len(d_2)+len(d_3)+len(d_4)+len(d_0)
        density.append(mixing_score)
print(len(cancer_cell))
print(len(immune_cell))
print(len(cancer_cell)/len(immune_cell))
print(np.mean(density))
