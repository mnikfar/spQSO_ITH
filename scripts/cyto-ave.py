import os
import pandas as pd
import numpy as np
import glob2 as glob
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
result_dir = './snapShots/'
patient_spatial = pd.read_csv(result_dir+'cyto_840.csv')
print(np.mean(patient_spatial ['IFNg']))
print(np.mean(patient_spatial ['IL_2']))