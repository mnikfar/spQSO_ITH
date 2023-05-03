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


data_folder = "graph2/" 
if not os.path.isdir(data_folder):
	os.makedirs(data_folder)
else:
	warnings.warn("Folder {} already exists, so it is removed first an then made".format(data_folder))
	shutil.rmtree(data_folder)
	os.makedirs(data_folder)
print('Saving graphs to '+ data_folder)
#data = pd.read_csv("solution.csv")
data_list=[]
color = ['red','black','blue','green']
pattern=['-','--','-.',':']
ylable1=["Stem",'Progenitor','Senesscent', 'Total']
ylable2=["Teff",'Tcyt','Texh','Treg-Default','Total']
legend=["I",'II','III', 'IV']
for k in range (4):
	data_list.append (pd.read_csv('stats_'+str(k)+".csv"))

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
num=0
#cancer cells
for i in range(5,9):
	num+=1
	figure(num=None, figsize=(width, length), dpi=DPI, facecolor='w', edgecolor='k')
	plt.ticklabel_format(axis='x',useMathText=True)
	#plt.ticklabel_format(axis='y',useMathText=True)
	#plt.ticklabel_format( axis='x', style='sci',scilimits=(0,0),useMathText=True)
	plt.ticklabel_format( axis='y', style='sci',scilimits=(0,0),useMathText=True)
	plt.tick_params(direction='in', length=6, width=2, colors='k',grid_color='k', grid_alpha=0.1,pad =Pad )
	plt.xlabel(data_list[0].columns[0]+' (day)',fontsize = 24, weight = 'bold',labelpad=10)
	plt.ylabel(ylable1[i-5] ,fontsize = 24 , weight = 'bold',labelpad=10)  
	plt.grid()
	for k in range (len(data_list)):
		if (i !=8):
			plt.plot(data_list[k][data_list[k].columns[0]]/4, data_list[k][data_list[k].columns[i]],pattern[k], c = color[k],label=legend[k],alpha=1.0,linewidth=LineWidth)
		else:
			cancer =data_list[k][data_list[k].columns[5]]+data_list[k][data_list[k].columns[6]]+data_list[k][data_list[k].columns[7]]
			plt.plot(data_list[k][data_list[k].columns[0]]/4, cancer,pattern[k], c = color[k],label=legend[k],alpha=1.0,linewidth=LineWidth)
	plt.legend(loc='best',fontsize = 15, ncol=1)
	plt.savefig(data_folder+str(num)+'.png', dpi=DPI)
	#plt.savefig(data_folder+data_list[k].columns[i]+'.svg',format='svg', dpi=DPI)
# (2) load this image into PIL
	if (tiff_save):
		png1 = Image.open(data_folder+data_list[k].columns[i]+'.png')
# (3) save as TIFF
		png1.save(data_folder+data_list[k].columns[i]+'.tiff')
	plt.close("all")
#T cells
for i in range(1,6):
	num+=1
	figure(num=None, figsize=(width, length), dpi=DPI, facecolor='w', edgecolor='k')
	plt.ticklabel_format(axis='x',useMathText=True)
	#plt.ticklabel_format(axis='y',useMathText=True)
	#plt.ticklabel_format( axis='x', style='sci',scilimits=(0,0),useMathText=True)
	plt.ticklabel_format( axis='y', style='sci',scilimits=(0,0),useMathText=True)
	plt.tick_params(direction='in', length=6, width=2, colors='k',grid_color='k', grid_alpha=0.1,pad =Pad )
	plt.xlabel(data_list[0].columns[0]+' (day)',fontsize = 24, weight = 'bold',labelpad=10)
	plt.ylabel(ylable2[i-1] ,fontsize = 24 , weight = 'bold',labelpad=10)  
	plt.grid()
	for k in range (len(data_list)):
		if (i !=5):
			plt.plot(data_list[k][data_list[k].columns[0]]/4, data_list[k][data_list[k].columns[i]],pattern[k], c = color[k],label=legend[k],alpha=1.0,linewidth=LineWidth)
		else:
			Tcell =data_list[k][data_list[k].columns[1]]+data_list[k][data_list[k].columns[2]]+data_list[k][data_list[k].columns[3]]+data_list[k][data_list[k].columns[4]]
			plt.plot(data_list[k][data_list[k].columns[0]]/4, Tcell,pattern[k], c = color[k],label=legend[k],alpha=1.0,linewidth=LineWidth)
	plt.legend(loc='best',fontsize = 15, ncol=1)
	plt.savefig(data_folder+str(num)+'.png', dpi=DPI)
	#plt.savefig(data_folder+data_list[k].columns[i]+'.svg',format='svg', dpi=DPI)
# (2) load this image into PIL
	if (tiff_save):
		png1 = Image.open(data_folder+data_list[k].columns[i]+'.png')
# (3) save as TIFF
		png1.save(data_folder+data_list[k].columns[i]+'.tiff')
	plt.close("all")
