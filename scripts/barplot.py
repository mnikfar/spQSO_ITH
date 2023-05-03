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
import glob2 as glob
from matplotlib.collections import LineCollection

data_folder = "bar/" 
if not os.path.isdir(data_folder):
    os.makedirs(data_folder)
else:
    warnings.warn("Folder {} already exists, so it is removed first an then made".format(data_folder))
    shutil.rmtree(data_folder)
    os.makedirs(data_folder)
print('Saving graphs to '+ data_folder)
#set Image Size, Resolustion, Fonts for plotting
length = 15 #inch
width = 8 #inch
DPI = 600
LabelSize = 20
axeslinewidth = 2.0
Pad = 10
LineWidth = 3.0
plotdata = pd.DataFrame({

    #"Case I":[0.0482,0.0351,0.0816],
    #"Case II":[0.2498,0.1853,0.3713],
    #"Case III":[0.0579,0.0414,0.0773],
    #"Case IV":[0.4382,0.3240,0.5422]  

    #"Case I":[0.0352,0.0213,0.0148],
    #"Case II":[0.0081,0.0046,0.0143],
    #"Case III":[0.0744,0.0443,0.0219],
    #"Case IV":[0.0578,0.0336,0.0884]  


    "M1: No-Treatment":[0.0482,0.2498,0.0579,0.4382],

    "M1: With With-Treatment":[0.0352,0.0081,0.0744,0.0578],

    "M2: No-Treatment":[0.0351,0.1853,0.0414,0.3240],

    "M2: With-Treatment":[0.0213,0.0081,0.0744,0.0578],

    "M3: No-Treatment":[0.0816,0.3713,0.0773,0.5422],

    "M3: With-Treatment":[0.0148,0.0143,0.0219,0.0884]    
    },

    index=["Case I", "Case II","Case III","Case IV"])
    #index=["M1", "M2","M3"])
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=LabelSize)
plt.rc('ytick', labelsize=LabelSize)
figure(num=None, figsize=(width, length), dpi=DPI, facecolor='w', edgecolor='k')
plotdata.plot(kind="bar",figsize=(length,width),color=['red', 'green','blue','yellow','pink','gray'],yerr=0.006,capsize= 2)
#plt.errorbar(p, w, yerr = wer,linestyle = 'None',c='green',capsize= 2)

plt.xticks(rotation=0, horizontalalignment="center",fontsize = 24, weight = 'bold')
plt.ylabel("Metrics",fontsize = 24, weight = 'bold',labelpad=10)
plt.grid()  
plt.legend(loc='best',fontsize = 15, ncol=1)
plt.savefig(data_folder+'bar.png', dpi=DPI)