import csv
import random
import sys
import os
import shutil
import argparse
import sys
import warnings
import pandas as pd
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from PIL import Image
from io import BytesIO
from scipy import special


data_folder = "sample/" 
if not os.path.isdir(data_folder):
	os.makedirs(data_folder)
else:
	warnings.warn("Folder {} already exists, so it is removed first an then made".format(data_folder))
	shutil.rmtree(data_folder)
	os.makedirs(data_folder)
print('Saving graphs to '+ data_folder)
#n = [0.0482,0.2498,0.0579,0.4382]
#w = [0.0352,0.0081,0.0744,0.0578]
#n= [0.0351,0.1853,0.0414,0.3240]
#w=[0.0213,0.0081,0.0744,0.0578]
 
n=[0.0816,0.3713,0.0773,0.5422]
ner= [number * 0.05 for number in n]
ner = 0.008


w=[0.0148,0.0143,0.0219,0.0884]
wer= [number * 0.05 for number in w]
wer = 0.008
num=4
p = range(1,num+1)
case = n+w
# sparate figure
#set Image Size, Resolustion, Fonts for plotting
width = 10 #inch
length = 8 #inch
DPI = 600
LabelSize = 20
axeslinewidth = 2.0
Pad = 10
LineWidth = 3.0
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=LabelSize)
plt.rc('ytick', labelsize=LabelSize)
plt.rcParams['axes.linewidth'] = axeslinewidth
figure(num=None, figsize=(width, length), dpi=DPI, facecolor='w', edgecolor='k')
plt.tick_params(direction='in', length=6, width=2, colors='k',grid_color='k', grid_alpha=0.1,pad =Pad )
plt.xlabel('Patients',fontsize = 24, weight = 'bold',labelpad=10)
plt.ylabel('M3' ,fontsize = 24 , weight = 'bold',labelpad=10)  
plt.xlim(0,num+1)
plt.grid()
plt.plot(p, n , c = 'red',label='No-Treatment',alpha=1.0,linestyle = 'None',marker='o')
plt.errorbar(p, n, yerr = ner,linestyle = 'None',c='red',capsize= 2)
plt.plot(p, w , c = 'green',label='With-Treatment',alpha=1.0,linestyle = 'None',marker='o')	
plt.errorbar(p, w, yerr = wer,linestyle = 'None',c='green',capsize= 2)
plt.legend(loc='best',fontsize = 15, ncol=1)
plt.savefig(data_folder+'test.png', dpi=DPI)
#plt.savefig(data_folder+data_list[k].columns[i]+'.svg',format='svg', dpi=DPI)
plt.close("all")


