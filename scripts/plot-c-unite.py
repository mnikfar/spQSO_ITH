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
data_folder = "graph/" 
if not os.path.isdir(data_folder):
	os.makedirs(data_folder)
else:
	warnings.warn("Folder {} already exists, so it is removed first an then made".format(data_folder))
	shutil.rmtree(data_folder)
	os.makedirs(data_folder)
print('Saving graphs to '+ data_folder)
data = pd.read_csv("stats_0.csv")
#set Image Size, Resolustion, Fonts for plotting
width = 12 #inch
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
legend1=["Stem",'Progenitor','Senesscent', 'Total']
legend2=["Teff",'Tcyt','Texh','Treg-Default','Total']
color = ['red','black','blue','green','brown']
pattern=['-','--','-.',':','-..']
#Cancer Cells
figure(num=None, figsize=(width, length), dpi=DPI, facecolor='w', edgecolor='k')
plt.ticklabel_format(axis='x',useMathText=True)
plt.ticklabel_format(axis='y',useMathText=True)
#plt.ticklabel_format( axis='x', style='sci',scilimits=(0,0),useMathText=True)
plt.ticklabel_format( axis='y', style='sci',scilimits=(0,0),useMathText=True)
plt.tick_params(direction='in', length=6, width=2, colors='k',grid_color='k', grid_alpha=0.1,pad =Pad )
plt.xlabel(data.columns[0]+' (day)',fontsize = 24, weight = 'bold',labelpad=10)  
plt.ylabel('Cell Number',fontsize = 24 , weight = 'bold',labelpad=10)
plt.grid()  
for i in range(5,9):
	if (i != 8):
		plt.plot(data[data.columns[0]]/4, data[data.columns[i]],pattern[i-5], c = color[i-5],label=legend1[i-5],alpha=1.0,linewidth=LineWidth)
	else:
		cancer =data[data.columns[5]]+data[data.columns[6]]+data[data.columns[7]]
		#cancer0 =data[data.columns[5]][0]+data[data.columns[6]][0]+data[data.columns[7]][0]
		#cancer = (cancer1-cancer0)/cancer0
		plt.plot(data[data.columns[0]]/4, cancer,pattern[i-5], c = color[i-5],label=legend1[i-5],alpha=1.0,linewidth=LineWidth)
	plt.legend(loc='best',fontsize = 15, ncol=1)
plt.savefig(data_folder+'C.png', dpi=DPI)
#plt.savefig(data_folder+data.columns[i]+'.svg',format='svg', dpi=DPI)
# (2) load this image into PIL
if (tiff_save):
	png1 = Image.open(data_folder+data.columns[i]+'.png')
# (3) save as TIFF
	png1.save(data_folder+data.columns[i]+'.tiff')
	png1.close
plt.close("all")

#T Cells
figure(num=None, figsize=(width, length), dpi=DPI, facecolor='w', edgecolor='k')
plt.ticklabel_format(axis='x',useMathText=True)
#plt.ticklabel_format(axis='y',useMathText=True)
#plt.ticklabel_format( axis='x', style='sci',scilimits=(0,0),useMathText=True)
plt.ticklabel_format( axis='y', style='sci',scilimits=(0,0),useMathText=True)
plt.tick_params(direction='in', length=6, width=2, colors='k',grid_color='k', grid_alpha=0.1,pad =Pad )
plt.xlabel(data.columns[0]+' (day)',fontsize = 24, weight = 'bold',labelpad=10)  
plt.ylabel('Cell Number',fontsize = 24 , weight = 'bold',labelpad=10)
plt.grid()  
for i in range(1,6):
	if (i != 5):
		plt.plot(data[data.columns[0]]/4, data[data.columns[i]],pattern[i-1], c = color[i-1],label=legend2[i-1],alpha=1.0,linewidth=LineWidth)
	else:
		Tcell =data[data.columns[1]]+data[data.columns[2]]+data[data.columns[3]]+data[data.columns[4]]
		#cancer0 =data[data.columns[5]][0]+data[data.columns[6]][0]+data[data.columns[7]][0]
		#cancer = (cancer1-cancer0)/cancer0
		plt.plot(data[data.columns[0]]/4, Tcell,pattern[i-1], c = color[i-1],label=legend2[i-1],alpha=1.0,linewidth=LineWidth)
	plt.legend(loc='best',fontsize = 15, ncol=1)
plt.savefig(data_folder+'T.png', dpi=DPI)
#plt.savefig(data_folder+data.columns[i]+'.svg',format='svg', dpi=DPI)
# (2) load this image into PIL
if (tiff_save):
	png1 = Image.open(data_folder+data.columns[i]+'.png')
# (3) save as TIFF
	png1.save(data_folder+data.columns[i]+'.tiff')
	png1.close
plt.close("all")

