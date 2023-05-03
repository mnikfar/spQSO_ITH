import numpy
import os
import pandas as pd
import geopandas
import pysal
import glob2 as glob
import numpy as np
import seaborn
import contextily
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from matplotlib.collections import LineCollection
from pointpats import centrography
from pointpats import (distance_statistics,QStatistic,random,PointPattern)
x0=242
z0=242
l=50
result_dir = './snapShots/'
for k in range (60,61,60):
    patient_spatial = pd.read_csv(result_dir+'score_{}.csv'.format(k))
    l_0 = []
    l_1 = []
    l_2 = []
    for y in range(5,6,1): 
        patient_slice = patient_spatial.loc[patient_spatial['y'] == y]
        cancer_cell = patient_slice.loc[patient_slice['Type'] == 1]
        n=0
        for i in range(len(cancer_cell)):
            x = cancer_cell.iloc[i]['x']
            z = cancer_cell.iloc[i]['z']
            if (x > x0 and x < x0+l and z > z0 and z < z0+l):
                n+=1
                cell_0 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z)]
                cell_1 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z+1)]
                cell_2 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z-1)]
                cell_3 = patient_slice.loc[(patient_slice['x'] == x+1) & (patient_slice['z'] == z)]
                cell_4 = patient_slice.loc[(patient_slice['x'] == x+1) & (patient_slice['z'] == z-1)]
                cell_5 = patient_slice.loc[(patient_slice['x'] == x+1) & (patient_slice['z'] == z+1)]
                cell_6 = patient_slice.loc[(patient_slice['x'] == x-1) & (patient_slice['z'] == z)]
                cell_7 = patient_slice.loc[(patient_slice['x'] == x-1) & (patient_slice['z'] == z-1)]
                cell_8 = patient_slice.loc[(patient_slice['x'] == x-1) & (patient_slice['z'] == z+1)]
                
                cell_9 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z+2)]
                cell_10 = patient_slice.loc[(patient_slice['x'] == x) & (patient_slice['z'] == z-2)]
                cell_11 = patient_slice.loc[(patient_slice['x'] == x+2) & (patient_slice['z'] == z)]
                cell_12 = patient_slice.loc[(patient_slice['x'] == x+2) & (patient_slice['z'] == z-2)]
                cell_13 = patient_slice.loc[(patient_slice['x'] == x+2) & (patient_slice['z'] == z+2)]
                cell_14 = patient_slice.loc[(patient_slice['x'] == x-2) & (patient_slice['z'] == z)]
                cell_15 = patient_slice.loc[(patient_slice['x'] == x-2) & (patient_slice['z'] == z-2)]
                cell_16 = patient_slice.loc[(patient_slice['x'] == x-2) & (patient_slice['z'] == z+2)]   
            
                d_0= cell_0.loc[cell_0['Type'] != 1]
                d_1= cell_1.loc[cell_1['Type'] != 1]
                d_2= cell_2.loc[cell_2['Type'] != 1]
                d_3= cell_3.loc[cell_3['Type'] != 1]
                d_4= cell_4.loc[cell_4['Type'] != 1]
                d_5= cell_5.loc[cell_5['Type'] != 1]
                d_6= cell_6.loc[cell_6['Type'] != 1]
                d_7= cell_7.loc[cell_7['Type'] != 1]
                d_8= cell_8.loc[cell_8['Type'] != 1]
                
                d_9= cell_9.loc[cell_9['Type'] != 1]
                d_10= cell_10.loc[cell_10['Type'] != 1]
                d_11= cell_11.loc[cell_11['Type'] != 1]
                d_12= cell_12.loc[cell_12['Type'] != 1]
                d_13= cell_13.loc[cell_13['Type'] != 1]
                d_14= cell_14.loc[cell_14['Type'] != 1]                 
                d_15= cell_15.loc[cell_15['Type'] != 1] 
                d_16= cell_16.loc[cell_16['Type'] != 1]  
    
                n_0=len(d_0)
                n_1=len(d_1)+len(d_2)+len(d_3)+len(d_4)+len(d_5)+len(d_6)+len(d_7)+len(d_8)
                n_2=len(d_9)+len(d_10)+len(d_11)+len(d_12)+len(d_13)+len(d_14)+len(d_15)+len(d_16)
                
                n_0=n_0
                n_1=n_0+n_1
                n_2=n_1+n_2
                
                l_0.append(n_0)
                l_1.append(n_1)
                l_2.append(n_2)
    L_0=np.mean(l_0)
    L_1=np.mean(l_1)
    L_2=np.mean(l_2)
R = np.arange(0,50.5,0.5)
L = []
for r in R:
    if (r == 0):
        L.append(0)
    if (r <= 10 and r !=0):
        L.append(1-np.exp(-L_0))
    if (r > 10 and r <= 30):
        L.append(1-np.exp(-L_1))
    if (r > 30 and r <= 50):
        L.append(1-np.exp(-L_2))
          
plt.plot(R,L)
plt.show()
AUC = 10*1-np.exp(-L_0)+20*(1-np.exp(-L_1))+20*(1-np.exp(-L_2))
print(AUC)