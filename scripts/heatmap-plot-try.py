#Required packages
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
import seaborn as sns

#set Image Size, Resolustion, Fonts for plotting
width = 10 #inch
length = 8 #inch
DPI = 600
tiff_save = False
LabelSize = 4
axeslinewidth = 2.0
Pad = 12
LineWidth = 0.5
Fontsize = 12
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=LabelSize)
plt.rc('ytick', labelsize=LabelSize)
plt.rcParams['axes.linewidth'] = axeslinewidth
#plt.ticklabel_format( axis='y', style='sci',scilimits=(0,0),useMathText=True)


#*****************Saving and generating folders
data_folder = "heatmap/" 
if not os.path.isdir(data_folder):
	os.makedirs(data_folder)
else:
	warnings.warn("Folder {} already exists, so it is removed first an then made".format(data_folder))
	shutil.rmtree(data_folder)
	os.makedirs(data_folder)
print('Saving graphs to '+ data_folder)


#*****************Plot the data
file = "QSP_3.csv"
data = pd.read_csv(file)
data = data.drop(data.columns[[0, 7]], axis=1)
fig, ax = plt.subplots(figsize=(width,length))
sns.heatmap(data.corr(), center=0, cmap='Reds', vmin=-1, vmax=1, annot=False,xticklabels=1, yticklabels=1)
#ax.set_xticklabels(ax.get_xticklabels(), rotation=10)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.set_title('Correlation Matrix',fontdict={'fontsize':Fontsize}, pad=Pad)
plt.savefig(data_folder+'heat.png', dpi=DPI,bbox_inches='tight')
plt.savefig(data_folder+'heat.svg',format='svg', dpi=DPI,bbox_inches='tight')
#plt.show()
plt.close("all")