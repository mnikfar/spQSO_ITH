#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 12:43:08 2020

Read simulation snapshots and create 2D drawings 
based on cell information

@author: changgong
"""


from PIL import Image, ImageDraw
from PIL import ImageOps


from pathlib import Path
import numpy as np

import lxml.etree as ET
import csv


#%%
# snapshot output enums
COL_X = 0
COL_Y = 1
COL_Z = 2
COL_CELLTYPE = 3
COL_CELLSTATE = 4

CELLTYEP_CANCER = 1
CELLTYEP_TCYT = 2
CELLTYEP_TREG = 3

STAINING_CD3 = 'CD3'
STAINING_CD8 = 'CD8'
STAINING_FoxP3 = 'FoxP3'
STAINING_PDL1 = 'PDL1'

############################### cross section #################################

# return an array of bool
# celldate: coordinates, type, state
# rad (r): radius of cells, unit: voxel
# loc: location of cut (integer, unit: voxel)
# dim: dimension of section (x,y,z = 0,1,2)
# thickness (h): thickness of cut, unit: voxel
# probability of cutting through cell can be calculated as:
#   p = min(2r+h, 1)
#   * in practice, no need to calculate min when 2r+h > 1
#   * if want to include all cells in on voxel layer, set h=1
def get_cross_section(celldata, loc, dim, thickness=1):
    selected = np.zeros(celldata.shape[0], dtype=bool)
    cell_in_layer = celldata[:,dim] == loc
    # 

    rad_layer = get_cell_radius(celldata[cell_in_layer, :])
    u01 = np.random.rand(rad_layer.shape[0])
    cell_in_slice = rad_layer*2+thickness > u01

    selected[cell_in_layer][cell_in_slice] = True
    selected[np.nonzero(cell_in_layer)[0][cell_in_slice]]=True
    return selected

# radius/voxel size
TYPE_RADIUS = {CELLTYEP_CANCER: .5,
       CELLTYEP_TCYT: .25,
       CELLTYEP_TREG: .25}

# return cell radius in an array, unit: voxel
def get_cell_radius(celldata):
    cell_type = celldata[:, COL_CELLTYPE]
    return np.array([TYPE_RADIUS[t] for t in cell_type])

#%%
############################### Basic shapes #################################

# n: number of edges; r_ab: long/short radius; 
# cent: center; rot: rotation, counterclock, degree
# return: list of (x, y)tuples
def get_polygon_ellipse(n=32, r_ab = [10,10], cent = [0, 0], rot = None):
    rad = np.array([i * 2 * np.pi / n for i in range(n)])
    pol = np.array([r_ab[0] * np.cos(rad), r_ab[1] * np.sin(rad)])
    return rotate_and_shift(pol, cent, rot)


def get_polygon_shape(pol0, r_ab = [10,10], cent = [0, 0], rot = None):
    pol = (pol0 * np.array(r_ab)).T
    return rotate_and_shift(pol, cent, rot)

def rotate_and_shift(pol, cent, rot):
     # rotate: non-traditional, because y is flipped in a image.
    if rot is not None:
        rr = rot/180*np.pi
        R = np.array([[np.cos(rr), np.sin(rr)],
                       [-np.sin(rr), np.cos(rr)]])
        pol = np.matmul(R, pol)
    # convert array to list of tuples    
    return (lambda pol:[(*p,) for p in pol])(pol.T + np.array(cent))
    
# load polygon vertices from svg file
def load_shapes(filename=None):
    if filename is None:
        filename = str(Path(__file__).parent/'cell_shapes.svg')
    tree = ET.parse(filename)
    root = tree.getroot()
    
    points_raw = [elem.attrib['points'] for elem in root if 'polygon' in elem.tag]
    polygon_arrays = [np.array(','.join(pol.split()).split(',')).astype('float').reshape(-1,2) for pol in points_raw]
    polygon_pool = []
    for i, pol in enumerate(polygon_arrays):
        xymax = np.max(pol, axis = 0)
        xymin = np.min(pol, axis = 0)
        polygon_pool += [((pol-xymin)/(xymax-xymin)-.5)*2]
    return polygon_pool

POLYGON_POOL = load_shapes()

#%%
############################## Cells as objects ##############################

#IHC_COLOR_BLUE = [200, 220, 250]
IHC_COLOR_BLUE = [150, 175, 250]
#COLOR_BROWN = [127, 100, 72]
IHC_COLOR_BROWN = [127, 80, 50]
IHC_COLOR_ALPHA_TISSUE = 20
IHC_COLOR_ALPHA_CELL = 50
IHC_COLOR_ALPHA_NUC = 80 

class cell:
    def __init__(self, rad, color):
        self.set_stain(color)
        self.r = rad
        return

    def set_stain(self, color):
        self._r = color[0]
        self._g = color[1]
        self._b = color[2]
        self._alpha = color[3]
        return

    # draw on to R/G/B channels
    def draw(self, imgs, cent, frac=1):
        r = self.r
        drawR = ImageDraw.Draw(imgs[0])
        drawG = ImageDraw.Draw(imgs[1])
        drawB = ImageDraw.Draw(imgs[2])
        drawA = ImageDraw.Draw(imgs[3])
        if self._r > 0:
            drawR.ellipse((cent[0]-r, cent[1]-r, cent[0]+r, cent[1]+r),
            fill = (self._r, 0, 0, 0))
        if self._g > 0:
            drawG.ellipse((cent[0]-r, cent[1]-r, cent[0]+r, cent[1]+r),
            fill = (0, self._g, 0, 0))
        if self._b > 0:
            drawB.ellipse((cent[0]-r, cent[1]-r, cent[0]+r, cent[1]+r),
            fill = (0, 0, self._b, 0))
        drawA.ellipse((cent[0]-r, cent[1]-r, cent[0]+r, cent[1]+r),
            fill = (0, 0, 0, int(self._alpha*frac)))
        return 

class cell_magnified(cell):

    def __init__(self, use_nuclei=True, use_env =True, use_shape=False):
        self.use_nuclei = use_nuclei
        self.use_env = use_env
        self.use_shape = use_shape
        return
    # f_deform: deformation factor, 1 is normal, 0 is no stretching, >1 is enhanced
    def set_shape(self, a, b, na, nb, mem_width = 1, f_deform = 1, f_env=1.2, n=32):
        self.n = n
        self.ra = a
        self.rb = b
        self.rna = na
        self.rnb = nb
        self.mem = mem_width
        self.f_deform = f_deform
        self.f_env = f_env
        return

class cell_magnified_ImF(cell_magnified):
    def __init__(self, use_nuclei=True, use_env =True, use_shape=False):
        super().__init__(use_nuclei, use_env, use_shape)
        self.set_stain()

    def set_stain(self, nuc = None, cyt = None, mem = None):
       
        self.stain_nuc = self.stain_cyt = self.stain_mem = False
        self.stain_nuc_rgba = self.stain_cyt_rgba = self.stain_mem_rgba = [0,0,0,0] 
        if nuc is not None:
            self.stain_nuc = True
            self.stain_nuc_rgba = nuc 
        if cyt is not None:
            self.stain_cyt = True
            self.stain_cyt_rgba = cyt
        if mem is not None:
            self.stain_mem = True
            self.stain_mem_rgba = mem
        return

    def draw(self, imgs, cent, ncent, fs=np.array([1,1]), rot = 0, f_stain = 1):
        
        r_ab = np.array([self.ra, self.rb])*np.exp(np.log(fs)*self.f_deform)
        r_nab = np.array([self.rna, self.rnb])*np.exp(np.log(fs)*self.f_deform)

        drawR = ImageDraw.Draw(imgs[0])
        drawG = ImageDraw.Draw(imgs[1])
        drawB = ImageDraw.Draw(imgs[2])
        
        # cell and  background tissue: 
        if self.use_shape:
            pol = POLYGON_POOL[int(rot)%3]
            
        if self.use_shape:
            pol_cell = get_polygon_shape(pol, r_ab = r_ab, cent = cent, rot = rot)
        else:
            pol_cell = get_polygon_ellipse(self.n, r_ab = r_ab, cent = cent, rot = rot)

        c = self.stain_cyt_rgba
        c[0] !=0 and drawR.polygon(pol_cell, fill =(c[0], 0, 0, c[3]))
        c[1] !=0 and drawG.polygon(pol_cell, fill =(0, c[1], 0, c[3]))
        c[2] !=0 and drawB.polygon(pol_cell, fill =( 0, 0, c[2],c[3]))

        # draw nuclei
        if self.use_nuclei:
            pol_nuc = get_polygon_ellipse(self.n, r_ab = r_nab, cent = ncent, rot = rot)
            c = self.stain_nuc_rgba
            c[0] !=0 and drawR.polygon(pol_nuc, fill =(c[0], 0, 0, c[3]))
            c[1] !=0 and drawG.polygon(pol_nuc, fill =(0, c[1], 0, c[3]))
            c[2] !=0 and drawB.polygon(pol_nuc, fill =(0, 0, c[2], c[3]))
        
        # draw membrane
        if self.stain_mem:
            c = self.stain_mem_rgba
            c[0] !=0 and drawR.line([*pol_cell, pol_cell[0]], fill = (c[0], 0, 0, c[3]), width = self.mem)
            c[1] !=0 and drawG.line([*pol_cell, pol_cell[0]], fill = (0, c[1], 0, c[3]), width = self.mem)
            c[2] !=0 and drawB.line([*pol_cell, pol_cell[0]], fill = (0, 0, c[2], c[3]), width = self.mem)
        return

class cell_magnified_IHC(cell_magnified):
    def __init__(self, use_nuclei=True, use_env =True, use_shape=False):
        super().__init__(use_nuclei, use_env, use_shape)
        self.set_stain()
    
    def set_stain(self, nuc = None, cyt = None, mem = None, env = None):
        self.base_rgb = IHC_COLOR_BLUE
        self.stain_rgb = IHC_COLOR_BROWN
        
        self.nuc_base_a = IHC_COLOR_ALPHA_NUC 
        self.cyt_base_a = IHC_COLOR_ALPHA_CELL
        self.mem_base_a = IHC_COLOR_ALPHA_CELL
        self.env_base_a = IHC_COLOR_ALPHA_TISSUE
        
        self.nuc_a = IHC_COLOR_ALPHA_NUC 
        self.cyt_a = IHC_COLOR_ALPHA_CELL
        self.mem_a = IHC_COLOR_ALPHA_CELL
        self.env_a = IHC_COLOR_ALPHA_TISSUE
        
        self.nuc_k = 1 
        self.cyt_k = 1
        self.mem_k = 1
        self.env_k = 1
        
        self.stain_nuc = self.stain_cyt = self.stain_mem = self.stain_env = False
        self.stain_nuc_rgb = self.stain_cyt_rgb = self.stain_mem_rgb = self.stain_env_rgb = self.base_rgb
        
        if nuc is not None:
            self.stain_nuc = True
            self.stain_nuc_rgb = self.stain_rgb
            self.nuc_a = nuc[0]
            self.nuc_k = nuc[1]
        if cyt is not None:
            self.stain_cyt = True
            self.stain_cyt_rgb = self.stain_rgb
            self.cyt_a = cyt[0]
            self.cyt_k = cyt[1]
        if mem is not None:
            self.stain_mem = True
            self.stain_mem_rgb = self.stain_rgb
            self.mem_a = mem[0]
            self.mem_k = mem[1]
        if env is not None:
            self.stain_env = True
            self.stain_env_rgb = self.stain_rgb
            self.env_a = env[0]
            self.env_k = env[1]
        return
    # base_a [0, 255]: opacity of base color
    # stain_a: [0, 255]: opacity of staining
    def get_mix_stain(self, base_a, stain_a):
        alpha = 255 - ((255 - base_a) * (255 - stain_a) / 255)
        red   = (self.base_rgb[0] * (255 - stain_a) + self.stain_rgb[0] * stain_a) / 255
        green = (self.base_rgb[1] * (255 - stain_a) + self.stain_rgb[1] * stain_a) / 255
        blue  = (self.base_rgb[2] * (255 - stain_a) + self.stain_rgb[2] * stain_a) / 255
        return (int(red), int(green), int(blue), int(alpha))
    
    # cent: center of cell; ncent: center of nuclei; 
    # fs: factor for cell bound, rot: rotation, degree
    # f_env: darkness of environment staining
    # f_stain: level of staining
    def draw(self, img, cent, ncent, fs=np.array([1,1]), rot = 0, f_stain = 1):
        
        _trans = Image.new('RGBA', img.size, (0,0,0,0))
        _drw = ImageDraw.Draw(_trans, 'RGBA')
        
        r_ab = np.array([self.ra, self.rb])*np.exp(np.log(fs)*self.f_deform)
        r_nab = np.array([self.rna, self.rnb])*np.exp(np.log(fs)*self.f_deform)
        
        
        # cell and  background tissue: 
        if self.use_shape:
            pol = POLYGON_POOL[int(rot)%3]
            
        if self.use_env:
            if self.use_shape:
                pol_env = get_polygon_shape(pol, r_ab = [self.f_env * x for x in r_ab], cent = cent, rot = rot)
            else:
                pol_env = get_polygon_ellipse(self.n, r_ab = [self.f_env * x for x in r_ab] , cent = cent, rot = rot)
            #_drw.polygon(pol_env, fill = (*self.stain_env_rgb,  int(f_stain* self.env_a)))
            _drw.polygon(pol_env, fill = self.get_mix_stain(self.env_base_a, self.env_a * f_stain**self.env_k))
                
        if self.use_shape:
            pol_cell = get_polygon_shape(pol, r_ab = r_ab, cent = cent, rot = rot)
        else:
            pol_cell = get_polygon_ellipse(self.n, r_ab = r_ab, cent = cent, rot = rot)
        #_drw.polygon(pol_cell, fill = (*self.stain_cyt_rgb,  self.cyt_a))
        _drw.polygon(pol_cell, fill = self.get_mix_stain(self.cyt_base_a, self.cyt_a*f_stain**self.cyt_k))

        # draw nuclei
        if self.use_nuclei:
            pol_nuc = get_polygon_ellipse(self.n, r_ab = r_nab, cent = ncent, rot = rot)
            #_drw.polygon(pol_nuc, fill = (*self.stain_nuc_rgb,  self.nuc_a) )
            _drw.polygon(pol_nuc, fill = self.get_mix_stain(self.nuc_base_a, self.nuc_a*f_stain**self.nuc_k))
        
        # draw membrane
        if self.stain_mem:
            _drw.line([*pol_cell, pol_cell[0]], 
                      #fill =  (*self.stain_mem_rgb,  int(f_stain*2/(1+f_stain) * self.mem_a)),
                      fill =  self.get_mix_stain(self.mem_base_a, self.mem_a*f_stain**self.mem_k),
                      width = self.mem)
                
            
        img.paste(Image.alpha_composite(img, _trans))
        
        return
    
    
##################################### Draw slide: ImF ###################################
## zoom out
def create_slide_ImF(cells, dim, res = 1):
    cancer_dim = 10 * res;
    t_dim = 8 * res;

    CancerCell = cell(cancer_dim, [0, 0, 255, 255])
    Tcyt = cell(t_dim, [0, 255, 0, 255])
    Treg = cell(t_dim, [255, 0, 0, 255])

    imgR = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,0))
    imgG = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,0))
    imgB = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,0))
    imgA = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,255))

    cell_proto = {
        CELLTYEP_CANCER: CancerCell,
        CELLTYEP_TCYT: Tcyt, 
        CELLTYEP_TREG: Treg }
    
    num_cell= len(cells)
    r_c = np.random.normal(scale = 3, size=(num_cell, 2))
    
    for i, c in enumerate(cells):
        cell_proto[c[2]].draw((imgR, imgG, imgB, imgA), (c[[0, 1]] + r_c[i])*res)

    rgba = np.array(imgR)+np.array(imgG)+np.array(imgB)+np.array(imgA)
    img = Image.fromarray(np.uint8(rgba))

    img = ImageOps.flip(img)     
    return img

## zoom in
def create_slide_ImF_closeup(cells, dim, res = 1, use_shape=False):
    
    # cells
    cancer_dim = np.array([12, 10, 5, 5])*res
    t_dim = np.array([5, 6, 3, 3])*res
     
    CancerCell = cell_magnified_ImF(use_nuclei = True, use_env = True, use_shape = use_shape)
    CancerCell.set_shape(*cancer_dim, mem_width = int(3*res), f_deform=1, f_env =1.5, n=16)
    CancerCell.set_stain(nuc = [0,0,255,255], cyt = [0,0,175,255], mem = None)
    
    Tcyt = cell_magnified_ImF(use_nuclei = True)
    Tcyt.set_shape(*t_dim, 0, f_deform=.2, n=16)
    Tcyt.set_stain(cyt = [0,200,0,255])
    
    Treg= cell_magnified_ImF(use_nuclei = True)
    Treg.set_shape(*t_dim, 0, f_deform=.2, n=16)
    Treg.set_stain(cyt = [200,0,0,255])

    imgR = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,255))
    imgG = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,255))
    imgB = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,255))

    cell_proto = {
        CELLTYEP_CANCER: CancerCell,
        CELLTYEP_TCYT: Tcyt, 
        CELLTYEP_TREG: Treg }

    num_cell= len(cells)
    r_c = np.random.normal(scale = 4, size=(num_cell, 2))
    r_nc = np.random.normal(scale = 1, size=(num_cell, 2))
    r_s = np.random.normal(scale = .25, size=(num_cell, 2))
    r_s2 = np.clip(np.exp(r_s), a_min=.5,a_max=1.5)
    r_r = np.random.rand(num_cell)
    
    for i, c in enumerate(cells):
        cent = np.array(c[[0, 1]] + r_c[i])*res
        ncent = cent + r_nc[i]*res
        stretch = r_s2[i]
        frac = c[3]
        rot = r_r[i] * 360
        cell_proto[c[2]].draw((imgR, imgG, imgB), cent, ncent, stretch, rot=rot, f_stain = frac)

    rgba = np.array(imgR)+np.array(imgG)+np.array(imgB)
    img = Image.fromarray(np.uint8(rgba))

    img = ImageOps.flip(img)     
    return img

##################################### Draw slide: IHC ###################################
## zoom out
def create_slide_IHC(cells, dim, res = 1, staining=None):

    color_no_stain = [*IHC_COLOR_BLUE, 127]        
    color_stain = [*IHC_COLOR_BROWN, 255]        
    cancer_dim = 10 * res;
    t_dim = 5 * res;
    CancerCell = cell(cancer_dim, color_no_stain)
    CancerCell0 = cell(cancer_dim, color_no_stain)
    Tcyt = cell(t_dim, color_no_stain)
    Treg = cell(t_dim, color_no_stain)
    T0 = cell(t_dim, color_no_stain)

    num_cell= len(cells)
    r_c = np.random.normal(scale = 3, size=(num_cell, 2))
    frac = np.ones(num_cell)
    khill = .5
    nhill = 4
    th_pdl1 = 0.2
        
    if staining == STAINING_CD3:
        Tcyt.set_stain(color_stain)
        Treg.set_stain(color_stain)
    elif staining == STAINING_CD8:
        Tcyt.set_stain(color_stain)
    elif staining == STAINING_FoxP3:
        Treg.set_stain(color_stain)
    elif staining == STAINING_PDL1:
        CancerCell.set_stain(color_stain)
        Tcyt.set_stain(color_stain)
        Treg.set_stain(color_stain)
        frac = cells[:,3]**nhill/(khill**nhill+cells[:,3]**nhill)

    imgR = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,0))
    imgG = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,0))
    imgB = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,0))
    imgA = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (0,0,0,0))

    cell_proto = {
        CELLTYEP_CANCER: CancerCell,
        CELLTYEP_TCYT: Tcyt, 
        CELLTYEP_TREG: Treg }
    cell_bg = {
        CELLTYEP_CANCER: CancerCell0,
        CELLTYEP_TCYT: T0, 
        CELLTYEP_TREG: T0
        }

    #fo = open('IHC.csv', 'w')
    for i, c in enumerate(cells):
        #fo.write('{},{},{},{}\n'.format(c[3], f, int(127*c[3]), int(127*f)))
        if frac[i]> th_pdl1:
            cell_proto[c[2]].draw((imgR, imgG, imgB, imgA), (c[[0, 1]] + r_c[i])*res, frac[i])
        else:
            cell_bg[c[2]].draw((imgR, imgG, imgB, imgA), (c[[0, 1]] + r_c[i])*res)

    #fo.close()

    rgba = np.array(imgR)+np.array(imgG)+np.array(imgB)+np.array(imgA)
    imga = Image.fromarray(np.uint8(rgba))

    img = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (255,255,255,255))
    img.paste(Image.alpha_composite(img, imga))

    img = ImageOps.flip(img)     
    return img

## zoom in
"""
cells:
    array: x, y, celltype, staining
dim: dimension (x, y) in microns
res: resolution, pixel per micron
staining: ['CD3', 'CD8', 'FoxP3', 'PDL1']
use_shape: use detail shape
"""
def create_slide_IHC_closeup(cells, dim, res = 1, staining=None, use_shape=False):
    
    khill = .5
    nhill = 4
    th_pdl1 = 0.2
    
    # cells
    cancer_dim = np.array([10, 8, 3, 3])*res
    t_dim = np.array([5, 4, 3, 3])*res
     
    CancerCell = cell_magnified_IHC(use_nuclei = True, use_env = True, use_shape = use_shape)
   
    CancerCell.set_shape(*cancer_dim, mem_width = int(3*res), f_deform=1, f_env =1.5, n=16)
    
    TCyt = cell_magnified_IHC(use_nuclei = True, use_env = True)
    TCyt.set_shape(*t_dim, 0, f_deform=.2, n=16)
    
    Treg= cell_magnified_IHC(use_nuclei = True, use_env = True)
    Treg.set_shape(*t_dim, 0, f_deform=.2, n=16)
    # staining
    if staining == STAINING_CD3:
        TCyt.set_stain(cyt=(220, 0), nuc=(220, 0))
        Treg.set_stain(cyt=(220, 0), nuc=(220, 0))
    elif staining == STAINING_CD8:
        TCyt.set_stain(cyt= (220, 0), nuc=(220, 0))
    elif staining == STAINING_FoxP3:
        Treg.set_stain(cyt=(220, 0), nuc=(220, 0))
    elif staining == STAINING_PDL1:
        CancerCell.set_stain(cyt = (100, 1), mem = (255, .5), env = (200, .5), nuc = (100, 1))
        TCyt.set_stain(cyt=(100, 0), nuc=(100, 0), env = (200, 1))
        Treg.set_stain(cyt=(100, 0), nuc=(100, 0), env = (200, 1))
    frac = cells[:,3]**nhill/(khill**nhill+cells[:,3]**nhill)
        
    img = Image.new('RGBA', (int(dim[0]*res), int(dim[1]*res)), (255,255,255,255))
    
    
    num_cell= len(cells)
    r_c = np.random.normal(scale = 4, size=(num_cell, 2))
    r_nc = np.random.normal(scale = 1, size=(num_cell, 2))
    r_s = np.random.normal(scale = .25, size=(num_cell, 2))
    r_s2 = np.clip(np.exp(r_s), a_min=.5,a_max=1.5)
    r_r = np.random.rand(num_cell)
    
    for i, c in enumerate(cells):
        cent = np.array(c[[0, 1]] + r_c[i])*res
        ncent = cent + r_nc[i]*res
        stretch = r_s2[i]
        rot = r_r[i] * 360
        
        cell_obj = None
        if c[2] == CELLTYEP_CANCER:
            cell_obj = CancerCell
        elif c[2] == CELLTYEP_TCYT:
            cell_obj = TCyt
        elif c[2] == CELLTYEP_TREG:
            cell_obj = Treg
        else :
            pass
        cell_obj.draw(img, cent, ncent, stretch, rot=rot, f_stain = frac[i])
    
    # flip image along the y axis    
    img = ImageOps.flip(img)     
    return img

#%% Other utilities

# construct path to an snapshot output file
def get_path(working_dir, sim_dir, prefix, t):
    path = str(working_dir/sim_dir) + '/snapShots/{}{}.csv'.format(prefix, t)
    return path

# retrieve data from csv file
def get_csv_data(filename, header_line = 1):
    with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        # header
        header = ''
        for i in range(header_line):
            header = next(reader)
        # data
        data = np.asarray(list(reader))
        return header, data

# get cell data from cell output file
# crd: x, y, z, type, state, numeric
# extra: cell specific additional data, text
def getCellData(cell_file_path):
    #print(cellFileName)
    header, data = get_csv_data(cell_file_path, 1)
    crd = data[:,:5].astype(float)
    extra=np.array([x.replace('"','') for x in data[:,-1]])
    return crd, extra

# get PDL1/maxPDL1 from extra output
PDL1_syn_MAX = 54658#*1.3
get_PDL1 = np.vectorize(lambda x: float(x.split('|')[0])/PDL1_syn_MAX)

# mappting to IHC cell type from simulation cell type
cell_type_map = {0: None,
                 1: CELLTYEP_CANCER,
                 2: CELLTYEP_TCYT,
                 3: CELLTYEP_TREG
                 }

# convert to IHC cell type id from simulation output cell type; vectorized
def get_cell_type(type_id):
    return np.vectorize(cell_type_map.get)(type_id)

# staining intensity, Hill's equation
def get_intensity(kd, n, x):
    return 1 / ( 1 + (kd/x)**n )

# slightly randomize cell coordinates
def randomize_coord(cell, factor):
    coord_rand = cell[:,:3] + np.random.rand(cell.shape[0], 3) * factor 
    return coord_rand 

def get_sim_dir(grp, sample, treat, rep):
    return 'group_{}/sample_{}/treatment_{}/rep_{}'.format(grp, sample, treat, rep)

def get_staining_PDL1(celldata, cellextra):
    PDL1_syn = get_PDL1(cellextra)
    #PDL1_syn = PDL1_syn**4*500
    kd = .1
    #kd = np.median(PDL1_syn)
    n = 5
    staining = get_intensity(kd, n, PDL1_syn)
    #staining = PDL1_syn
    return staining

def get_staining_CD3(celldata, cellextra):
    staining = np.logical_or (get_cell_type(celldata[:,3])==CELLTYEP_TCYT, 
                      get_cell_type(celldata[:,3])==CELLTYEP_TREG)
    return staining

def get_staining_CD8(celldata, cellextra):
    staining = get_cell_type(celldata[:,3])==CELLTYEP_TCYT
    return staining

def get_staining_FoxP3(celldata, cellextra):
    staining = get_cell_type(celldata[:,3])==CELLTYEP_TREG
    return staining

get_staining = {STAINING_CD3: get_staining_CD3,
                STAINING_CD8: get_staining_CD8,
                STAINING_FoxP3: get_staining_FoxP3,
                STAINING_PDL1: get_staining_PDL1}

   
#%%
    
if (__name__ == '__main__'):

    cells = np.array([[30, 30, 1, .8], 
                      [40, 60, 2, .2],
                      [50, 50, 3, .5]])
    
    # distance ImF
    #img = create_slide_ImF(cells, [100,100], res = 5)
    # close up ImF 
    img = create_slide_ImF_closeup(cells, [100,100], res = 5, use_shape=True)

    # close up IHC
    # close up IHC
    #img = create_slide_IHC_closeup(cells, [100,100], res = 5, staining='PDL1', use_shape=True)

    img.show()