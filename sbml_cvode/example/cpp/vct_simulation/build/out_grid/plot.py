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


data_folder = "graph/" 
if not os.path.isdir(data_folder):
	os.makedirs(data_folder)
else:
	warnings.warn("Folder {} already exists, so it is removed first an then made".format(data_folder))
	shutil.rmtree(data_folder)
	os.makedirs(data_folder)
print('Saving graphs to '+ data_folder)
#data = pd.read_csv("solution.csv")
data_list=[]
color = []
for k in range (15):
	data_list.append (pd.read_csv(str(k)+".csv"))
	c = (random.random(), random.random(), random.random())
	color.append(c)



# sparate figure
#set Image Size, Resolustion, Fonts for plotting
width = 10 #inch
length = 8 #inch
DPI = 600
tiff_save = False
LabelSize = 20
axeslinewidth = 2.0
Pad = 10
LineWidth = 3.0
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=LabelSize)
plt.rc('ytick', labelsize=LabelSize)
plt.rcParams['axes.linewidth'] = axeslinewidth
for i in range(len(data_list[0].columns)):
	if (i!=0):
		figure(num=None, figsize=(width, length), dpi=DPI, facecolor='w', edgecolor='k')
		for k in range (len(data_list)):
			plt.ticklabel_format(axis='x',useMathText=True)
			#plt.ticklabel_format(axis='y',useMathText=True)
			#plt.ticklabel_format( axis='x', style='sci',scilimits=(0,0),useMathText=True)
			plt.ticklabel_format( axis='y', style='sci',scilimits=(0,0),useMathText=True)
			plt.tick_params(direction='in', length=6, width=2, colors='k',grid_color='k', grid_alpha=0.1,pad =Pad )
			plt.xlabel(data_list[k].columns[0]+' (month)',fontsize = 24, weight = 'bold',labelpad=10)  
			plt.ylabel(data_list[k].columns[i] ,fontsize = 24 , weight = 'bold',labelpad=10)
			plt.grid()
			plt.plot(data_list[k][data_list[k].columns[0]]/86400/30, data_list[k][data_list[k].columns[i]],'-', c = color[k],label=k+1,alpha=1.0,linewidth=LineWidth)
			plt.legend(loc='best',ncol=3)
			#plt.plot(data[data.columns[0]]/86400/30, data[data.columns[3]],'-', c = color,label='',alpha=1.0,linewidth=LineWidth)
			plt.savefig(data_folder+data_list[k].columns[i]+'.png', dpi=DPI)
			#plt.savefig(data_folder+data_list[k].columns[i]+'.svg',format='svg', dpi=DPI)
# (2) load this image into PIL
			if (tiff_save):
				png1 = Image.open(data_folder+data_list[k].columns[i]+'.png')
# (3) save as TIFF
				png1.save(data_folder+data_list[k].columns[i]+'.tiff')
		plt.close("all")