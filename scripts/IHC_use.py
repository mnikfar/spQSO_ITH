"""
Created on Thu May 14 13:58:58 2020

@author: Chang
"""

import sys
import numpy as np
from pathlib import Path
import importlib
import spQSP_histpath_visual as ihc
importlib.reload(ihc)


#%%
def draw_snapshot(input_filename, output_dir, serial):
    #%%
    out_imf= output_dir/'imf'
    out_imf_large = output_dir/'imf_zoom'

    out_imf.mkdir(parents=True, exist_ok=True)
    out_imf_large.mkdir(parents=True, exist_ok=True)

    #filename = get_path(working_dir, sim_dir, 'cell_', t_slice)
    celldata, cellextra = ihc.getCellData(input_filename)

    #%% slice
    selected = ihc.get_cross_section(celldata, 5, 1, 1)

    cells = np.zeros((np.count_nonzero(selected), 4))
    cells[:,:2] = celldata[selected][:,[0,2]]* vox
    cells[:,2] = celldata[selected,3]
    #%% imF full
    img = ihc.create_slide_ImF(cells, domain_size, res = res_wsi)
    img.convert('RGB').save(str(out_imf)\
     + '/imf_wsi_{}.png'.format(serial), "PNG")
    #% zoom in
    x1,z1 = x0+x_length,z0+z_length
    x_cut = np.logical_and(cells[:,0]<x1, cells[:,0]>x0) 
    z_cut = np.logical_and(cells[:,1]<z1, cells[:,1]>z0) 
    mag_window = np.logical_and(x_cut, z_cut)
    #%
    #% ImF_closeup
    img = ihc.create_slide_ImF_closeup(cells[mag_window]-np.array([x0,z0,0,0]).T, 
                                   [x_length, z_length], res = res_zoom, use_shape=True)
    img.convert('RGB').save(str(out_imf_large)\
     + '/imf_{}.png'.format(serial), "PNG")

    return


#%%
# input
working_dir = Path('.')
# params
vox = 20 # voxel size, micron
# wsi
domain_size = [10000, 10000] # micron
res_wsi = 0.2 # pixel per micron
res_zoom = 1# pixel per micron
x0,z0,= 2800, 2800# starting corner, microns 
x_length, z_length = 4500, 4500 # zoom in window, microns 
filename = 'score_360.csv'
print('Drawing for file: {}'.format(filename))
output_dir = working_dir/'sample_240'
draw_snapshot(filename, output_dir, 'cell_240')