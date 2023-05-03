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
from io import BytesIO
from scipy import special

#**************************Saving and generating folders*********************
data_folder = "Tecplot/" 
if not os.path.isdir(data_folder):
	os.makedirs(data_folder)
else:
	warnings.warn("Folder {} already exists, so it is removed first an then made".format(data_folder))
	shutil.rmtree(data_folder)
	os.makedirs(data_folder)
print('Saving data to '+ data_folder)
#***************************Convert data***************************************
num_index = 100
step = 1
SizeX = 100
SizeY = 100
SizeZ = 100
nc=SizeX*SizeY*SizeZ 
f1 = open(data_folder+"grid.dat", "w")
f1.write("TITLE = \"Scalar Field\"\n")
f1.write("VARIABLES = \"x\", \"y\", \"z\"\n")
f1.write("Zone T = \"Frame 0\", I = {}, J = {}, K = {}\n".format(SizeX+1,SizeY+1,SizeZ+1))
for k in range(SizeZ+1):
	for j in range(SizeY+1):
		for i in range(SizeX+1):
				f1.write("{} {} {}\n".format(i,j,k))
f1.close()
for kk in range (0, num_index+1, step):
	num = 0
	core= pd.read_csv("snapShots/"+"cyto_"+str(kk)+".csv")
	f1 = open(data_folder+"cyto_"+str(kk)+".dat", "w")
	f1.write("TITLE = \"Scalar Field\"\n")
	f1.write("VARIABLES = \"x\", \"y\", \"z\", \"IFNg\", \"IL2\"\n")
	f1.write("Zone T = \"Frame 0\", I = {}, J = {}, K = {}\n".format(SizeX,SizeY,SizeZ))
	for k in range(SizeZ):
		for j in range(SizeY):
			for i in range(SizeX):
				f1.write("{} {} {} {} {}\n".format(i+0.5,j+0.5,k+0.5,core.iloc[num,0],core.iloc[num,1]))
				#f2.write("{} {} {} {} {}\n".format(i,j,k,margin.iloc[num,0],margin.iloc[num,1]))
				num=num+1
	f1.close()



