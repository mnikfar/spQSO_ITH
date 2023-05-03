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


data_folder = "Metrics/" 
if not os.path.isdir(data_folder):
	os.makedirs(data_folder)
else:
	warnings.warn("Folder {} already exists, so it is removed first an then made".format(data_folder))
	shutil.rmtree(data_folder)
	os.makedirs(data_folder)
print('Saving graphs to '+ data_folder)

Time = [15,30,45,60,75,90,105,120,135,150,165,180] 
M1n = [1,1,1,1,1,1,1,1,1,1,1,1]
M1ner= [number * 0.00 for number in M1n]


#M2n = [0.820783216783217,0.8181818181818182,0.8239209388742099,0.8232742257742258,0.8273643023643026,0.8290622421057204,0.8272668507962625,0.8084407938999778,0.8156317134193242,0.824482096850518,0.8359669991025923,0.8039363861944506]
#M2ner= [number * 0.005 for number in M2n]


#M3n = [0.9388642389031939,0.9434511710149045,0.9585937448656856,0.9393820678814518,0.95483398694108,0.9366487265257744,0.9683984469660692,0.9812422584243229,0.9755091513469198,0.9536251076179925,0.9690419761526357,0.9519217266471056]
#M3ner= [number * 0.005 for number in M3n]

M1w = [1,1,1,1,1,1,1,0.9930394431554525,0.9615384615384616,0.9393139841688655,0.9434523809523809,0.98]
M1wer= [number * 0.004 for number in M1w]

#M2w = [0.820783216783217,0.8181818181818182,0.8239209388742099,0.8232742257742258,0.7488640080157938,0.7053725758787771,0.69112959604948,0.6739303289901221,0.6575891918538977,0.6353290010642952,0.6344217378497158,0.654608795496646]
#M2wer= [number * 0.005 for number in M2w]

#M3w = [0.9388642389031939,0.9434511710149045,0.9585937448656856,0.9393820678814518,0.8646788535953924,0.7744853503994168,0.7482002694195887,0.7683833757481402,0.7570750800802963,0.761587888702739,0.7454939950852266,0.7926859081261712]
#M3wer= [number * 0.005 for number in M3w]




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
plt.xlabel('Time (day)',fontsize = 24, weight = 'bold',labelpad=10)
plt.ylabel('Mixing Score' ,fontsize = 24 , weight = 'bold',labelpad=10) 
plt.xlim(15,185)
plt.ylim(0.7,1)
plt.xticks(np.arange(15, 185, 15.0))
plt.yticks(np.arange(0.6, 1.1, 0.1))
plt.grid()
plt.plot(Time, M1n , c = 'red',label='No Treatment',alpha=1.0,linestyle = '-',marker='o')
plt.errorbar(Time, M1n, yerr = M1ner,linestyle = 'None',c='red',capsize= 2)




plt.plot(Time, M1w , c = 'green',label='With Treatment',alpha=1.0,linestyle = '-',marker='s')	
plt.errorbar(Time, M1w, yerr = M1wer,linestyle = 'None',c='green',capsize= 2)



#plt.axvspan(1-0.25, 1+0.25, color='yellow', alpha=0.75, lw=0)
#plt.axvspan(3-0.25, 3+0.25, color='yellow', alpha=0.75, lw=0)
#plt.axvspan(5-0.25, 5+0.25, color='yellow', alpha=0.75, lw=0)
#plt.axvspan(7-0.25, 7+0.25, color='yellow', alpha=0.75, lw=0)
#plt.axvspan(8-0.25, 8+0.25, color='gray', alpha=0.75, lw=0)
#plt.axvspan(12-0.25, 12+0.25, color='pink', alpha=1, lw=0)



plt.legend(loc='best',fontsize = 15, ncol=1)
plt.savefig(data_folder+'test.png', dpi=DPI)
#plt.savefig(data_folder+data_list[k].columns[i]+'.svg',format='svg', dpi=DPI)
plt.close("all")


