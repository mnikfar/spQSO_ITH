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


data_folder = "graph3/" 
if not os.path.isdir(data_folder):
	os.makedirs(data_folder)
else:
	warnings.warn("Folder {} already exists, so it is removed first an then made".format(data_folder))
	shutil.rmtree(data_folder)
	os.makedirs(data_folder)
print('Saving graphs to '+ data_folder)
data = pd.read_csv("sim_results.csv")


case = 1 # 1 for subplots and 2 for separtes plots
if (case == 1):
# subplots
#set Image Size, Resolustion, Fonts for plotting 
	Nrows=3 
	Ncols=3
	width = 10 #inch
	length = 40 #inch
	DPI = 600
	tiff_save = False
	LabelSize = 4
	axeslinewidth = 1
	Pad = 1
	LineWidth = 1.0
	plt.rc('font', family='serif')
	plt.rc('xtick', labelsize=LabelSize)
	plt.rc('ytick', labelsize=LabelSize)
	plt.rcParams['axes.linewidth'] = axeslinewidth
	fig, ax = plt.subplots(nrows=Nrows, ncols=Ncols)
	i = 0
	for row in ax:
		i=i+1
		j = 0
		for col in row:
			j=j+1
			num=i*j
			if (num < len(data.columns)):
				r = random.random()
				b = random.random()
				g = random.random()
				color = (r, g, b)
				plt.ticklabel_format(axis='x',useMathText=True)
				col.ticklabel_format( axis='y', style='sci',scilimits=(0,0),useMathText=True)
				col.tick_params(direction='in', length=1, width=1, colors='k',grid_color='k', grid_alpha=0.1,pad =Pad)
				col.set_xlabel(data.columns[0]+' (month)' ,fontsize = LabelSize , weight = 'bold',labelpad=1)  
				col.set_ylabel(data.columns[num] ,fontsize = LabelSize , weight = 'bold',labelpad=1)
				col.grid()  
				col.plot(data[data.columns[0]]/86400/30, data[data.columns[num]],'-', c = color,label='',linewidth =axeslinewidth)
	plt.savefig(data_folder+'Results.png', dpi=DPI)
	plt.savefig(data_folder+'Results.svg',format='svg', dpi=DPI)
# (2) load this image into PIL
	if (tiff_save):
		png1 = Image.open('Results.png')
# (3) save as TIFF
		png1.save(data_folder+data.columns[i]+'Results.tiff')
		png1.close
	plt.close("all")
else:
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
	for i in range(len(data.columns)):
		if (i!=0):
			figure(num=None, figsize=(width, length), dpi=DPI, facecolor='w', edgecolor='k')
			r = random.random()
			b = random.random()
			g = random.random()
			color = (r, g, b)
			plt.ticklabel_format(axis='x',useMathText=True)
			#plt.ticklabel_format(axis='y',useMathText=True)
			#plt.ticklabel_format( axis='x', style='sci',scilimits=(0,0),useMathText=True)
			plt.ticklabel_format( axis='y', style='sci',scilimits=(0,0),useMathText=True)
			plt.tick_params(direction='in', length=6, width=2, colors='k',grid_color='k', grid_alpha=0.1,pad =Pad )
			plt.xlabel(data.columns[0]+' (month)',fontsize = 24, weight = 'bold',labelpad=10)  
			plt.ylabel(data.columns[i] ,fontsize = 24 , weight = 'bold',labelpad=10)
			plt.grid()  
			plt.plot(data[data.columns[0]]/86400/30, data[data.columns[i]],'-', c = color,label='',alpha=1.0,linewidth=LineWidth)
			plt.savefig(data_folder+data.columns[i]+'.png', dpi=DPI)
			#plt.savefig(data_folder+data.columns[i]+'.svg',format='svg', dpi=DPI)
# (2) load this image into PIL
			if (tiff_save):
				png1 = Image.open(data_folder+data.columns[i]+'.png')
# (3) save as TIFF
				png1.save(data_folder+data.columns[i]+'.tiff')
				png1.close
			plt.close("all")