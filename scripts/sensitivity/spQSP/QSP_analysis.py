# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 17:07:23 2020

Functions for model result analysis and visualization

@author: Chang

Index of functions:
    
* Data fetching and preprocessing
    * Read data from .csv files
    * Obtain LHS matrix
    * Preprocess raw output and obtain endpoints and other readouts
        * Tumor diameter
        * Tumor volume
        * Time to progression
        * Time to response
* Generating figures:
    * 2-way clustered heatmap (PRCC)
    * waterfall (Volume change, sorted)
    * percentage vs percentile (RECIST vs param/readout)
    * survival
    * spider
    * boxplot for group differences
    * ROC
    * Hazard ratio

"""
# Plotting
import seaborn as sns
# numerical
import csv 
# file handling
import numpy as np 
import pandas as pd
from itertools import groupby


#####################################################################
#
#   Data fetching and preprocessing
#
#####################################################################

"""
Read data from .csv files
Inputs: 
    filename: path to .csv file
    header line: number of rows for headers. The last one is outputed
    dtype: return data type, default to str
Outputs: 
    header: list of strings
    data: numpy 2D array
    
"""
def read_csv(filename, header_line = 1, dtype = str):
    with open(filename) as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        # header
        header = ''
        for i in range(header_line):
            header = next(reader)
        # data
        data = np.asarray(list(reader), dtype = dtype)
        return header, data    
"""
Obtain LHS matrix
Inputs:
    filename: path to .csv file. (first column: exp id; first row: Parameters)
Outputs:
    lhs: n x k, parameter values, n experiments vs k parameters
    exp_id: n, array of exp ids.
    param_name: k, list of parameters included
"""
def read_LHS(filename):
    header, data = read_csv(filename)
    lhs = data[:,1:]
    exp_id = data[:,0].astype(int)
    param_name = header[1:]
    return  lhs, exp_id, param_name


"""
Wrtite data to a csv file
"""
def write_csv(filename, data, header):
    with open(filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(header)
        csvwriter.writerows(data) 
    return

"""
Tumor diameter
Inputs: 
    num_T: CD8 + Treg
    num_C: Cancer cell number
    num_APC: APC number
    scalor: mulitplier to all cell numbers (for weighted QSP)
    f_vol_tum: interstitial fraction
"""
def get_tum_volume(num_T, num_C, num_APC, scalor=1, f_vol_tum = 0.718):
    vol_T = 3.1415926/6*(6.94/1000)**3
    vol_C = 3.1415926/6*(16.9/1000)**3
    vol_APC = 3.1415926/6*(9.14/1000)**3
    vol_tum = (vol_T * num_T + vol_C * num_C + vol_APC * num_APC) * scalor /(1-f_vol_tum)
    vol_tum = vol_tum.clip(min=0)
    return vol_tum

def get_tum_diameter(num_T, num_C, num_APC, scalor=1):
    vol_tum = get_tum_volume(num_T, num_C, num_APC, scalor)
    D_tum_app = (vol_tum/(np.pi/6))**(1/3)
    return D_tum_app
#####################################################################
#
#   Visualization
#
#####################################################################

"""
Seaborn clustermap plot. For PRCC
"""
#
def cluster_map(data, row_label, col_label, fig_size, annot = None,row_color_labels=None, col_color_labels=None,
                show_dendrogram = [True, True], **kwarg):
    
    df = pd.DataFrame(data=data, index = row_label, columns = col_label)
    g = sns.clustermap(df, annot = annot, fmt = '',
                       vmin=-1, vmax=1, cbar_kws={"ticks":[-1, -.5,  0, .5,  1]}, **kwarg)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, va = 'center')
    if row_color_labels is not None:
        row_colors = kwarg['row_colors']
        borders = np.cumsum([0] + [sum(1 for i in g) for k, g in groupby(row_colors)])
        for b0, b1, label in zip(borders[:-1], borders[1:], row_color_labels):
            g.ax_row_colors.text(-0.06, (b0 + b1) / 2, label, color='black', ha='right', va='center', rotation=90,
                                 transform=g.ax_row_colors.get_yaxis_transform())
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=-55, ha = 'left')
    if col_color_labels is not None:
        col_colors = kwarg['col_colors']
        borders = np.cumsum([0] + [sum(1 for i in g) for k, g in groupby(col_colors)])
        for b0, b1, label in zip(borders[:-1], borders[1:], col_color_labels):
            g.ax_col_colors.text((b0 + b1) / 2, 1.06, label, color='black', ha='center', va='bottom',
                                 transform=g.ax_col_colors.get_xaxis_transform())
    g.ax_row_dendrogram.set_visible(show_dendrogram[0])
    g.ax_col_dendrogram.set_visible(show_dendrogram[1])
    g.fig.set_size_inches(*fig_size)
    return g


"""
y: diameter change %
group id: each y belongs to one of k groups, with id of 0, 1, ..., k-1 
clist: color code for each group 
txt_x: location of RECIST annotation
"""
def waterfall(ax, y, group_id, clist, txt_x=[.05, .8], by_group=True, **kwarg):
    # sort
    num_case = len(y)
    rank = np.zeros(num_case, dtype=int)
    seq = np.arange(num_case)
    groups = np.unique(group_id)
    if by_group:
        ng = 0
        for g in groups:
            idx = (group_id == g)
            rank[ng:ng+np.sum(idx)] = seq[idx][np.flip(np.argsort(y[idx]))]
            ng += np.sum(idx)
    else:
        rank = seq[np.flip(np.argsort(y))]
    # barplot
    colors = np.asarray([clist[x] for x in group_id])
    ax.bar(range(num_case), y[rank]*100, color = colors[rank])
    # extra
    ax.hlines(-30, 0, num_case, linewidth=1, color='grey', ls=':')
    ax.hlines(20, 0, num_case, linewidth=1, color='grey', ls=':')
    ax.text(num_case*txt_x[0], -26, '-30 %', fontsize=12)
    ax.text(num_case*txt_x[1], 24, '20 %', fontsize=12)
    ax.set_xlim(0,num_case)
    ax.set_ylim(-100,100)
    ax.get_xaxis().set_ticks([])
    return


"""
alternative waterfall plot. show areas instead of bars
y: diameter change %
group id: each y belongs to one of k groups, with id of 0, 1, ..., k-1 
clist: color code for each group 
txt_x: location of RECIST annotation
"""
def waterfall_b(ax, y, group_id, clist, **kwarg):
    # sort
    num_case = len(y)
    rank = np.zeros(num_case, dtype=int)
    seq = np.arange(num_case)
    groups = np.unique(group_id)
    ng = 0
    for i, g in enumerate(groups):
        idx = (group_id == g)
        ng_new = ng+np.sum(idx)
        rank[ng:ng_new] = seq[idx][np.flip(np.argsort(y[idx]))]
        ax.fill_between(np.arange(ng, ng_new), 0, y[rank[ng:ng_new]]*100, facecolor=clist[i])
        ng = ng_new
    # extra
    ax.hlines(-30, 0, num_case, linewidth=.5, color='k', ls='--')
    ax.hlines(20, 0, num_case, linewidth=.5, color='k', ls='--')
    ax.set_xlim(0,num_case)
    ax.set_ylim(-100,100)
    ax.get_xaxis().set_ticks([])
    return

"""
plot fraction of different response vs quantile of examined variable
y: response [-1, +inf]
x: parameter/readout values
n_xbin: how many bins for x
seg_progressive: cutoff for levels of progression
seg_responsive: cutoff for levels of responsive ness
"""
def response_q(ax, y, x, n_xbin, 
               seg_progressive=[.2, .5, 1], seg_responsive=[-.3, -.9, -.99],
               color_R='tab:green', color_NR= 'tab:red', alpha = .15, **kwarg):
    # percentiles
    p = np.arange(0,n_xbin)*100/n_xbin
    print(p)
    p_mid = p + 50/n_xbin
    num_seg_progressive = len(seg_progressive)
    num_seg_responsive = len(seg_responsive)
    # 
    recist_pct = np.zeros((n_xbin, num_seg_progressive + num_seg_responsive))
    # RECIST fraction
    for i, p_low in enumerate(p):
        p_hi = p_low + 100/n_xbin
        x_min = np.percentile(x, p_low)
        x_max = np.percentile(x, p_hi)
        idx = np.logical_and(x >= x_min, x < x_max)
        for j, s in enumerate(seg_progressive):
            recist_pct[i, j] = sum(y[idx]<s)/sum(idx)*100
        for j, s in enumerate(seg_responsive):
            recist_pct[i, -j-1] = sum(y[idx]<s)/sum(idx)*100
    # line plot
    #ax.plot(p_mid, recist_pct[:,-1], linewidth=.5, color='tab:green')
    #ax.plot(p_mid, recist_pct[:,1], linewidth=.5, color='tab:red')
    ax.fill_between(p_mid, recist_pct[:,-1], recist_pct[:,0], color = 'white', alpha=alpha)
    for j, s in enumerate(seg_progressive):
        ax.fill_between(p_mid, recist_pct[:,j], 100, color = color_NR, alpha=alpha)
    for j, s in enumerate(seg_responsive):   
        ax.fill_between(p_mid, 0, recist_pct[:,-j-1], color = color_R, alpha=alpha)
    # extra
    ax.set_xlim(50/n_xbin,100-50/n_xbin)
    ax.set_ylim(0,100)
    return


"""
spider plot of tumor size change over time
ax: axis
t: time
y: volume change
colors: colors for line
"""
def spider(ax, t, y, colors, alpha, **kwarg):
    for i, c in enumerate(colors):
        yi = y[i]        
        line = ax.plot(t, yi, alpha=alpha, color = c)
    return

"""
survival plot
time_to_events: list by groups. Time when event occur for each individual. -1 indicate no event.
invert: True for percent event happened; False for remaining (survival, default)
"""
def survival(ax, t_max, time_to_events, invert=False, **kwarg):
    t_all = []
    pct_all = []
    inv = 1*invert
    for y in time_to_events:
        yt = y[y>0]
        n = len(yt)
        n_tot = len(y)
        yt.sort()
        _t = np.zeros(n*2+2)
        _y = np.zeros(n*2+2)
        _t[0], _y[0] = 0, 100
        for i, t in enumerate(yt):
            _t[i*2+1] = t
            _t[i*2+2] = t
            _y[i*2+1] = _y[i*2]
            _y[i*2+2] = (1-(i+1)/n_tot)*100
        _t[-1] = t_max
        _y[-1] = _y[-2]
        t_all += [_t]
        pct_all += [_y]
    for i, ti in enumerate(t_all):
        surv = 100*inv + (1-2*inv)*pct_all[i]
        ax.plot(ti, surv, **kwarg)
    ax.set_xlim([0, t_max])
    ax.set_ylim([-5, 105])
    return


"""
Forest plot
data: dataframe with columns: name, n, x, cl, cu
"""
import matplotlib.lines as mlines

def forest_plot(axs, data, xlabel, xlim, xbar, label, f_table=.5, fmt='g', lvl=2):
    n_row = len(data)
    Y = []
    X = []
    error_lower = []
    error_upper = []
    table_content = []
    for i, entry in data.iterrows():
        if np.isnan(entry['n']):
            # name only
            table_content.append(['$\\bf{}$'.format(entry['name'].split('.')[-1]).replace('_', '\_'), '', '', ''])
            pass
        else:
            y = n_row - i -.5
            Y.append(y)
            x = entry['x']
            X.append(x)
            cl = entry['cl']
            cu = entry['cu']
            error_lower.append(x-cl)
            error_upper.append(cu-x)
            if lvl == 1:
                row_name = entry['name'].split('.')[-1]
            else:
                row_name = entry['name']
            #print(Y, X, error_lower, error_upper)
            table_content.append([row_name, int(entry['n']), 
                                 ('{:'+fmt+'}').format(x), 
                                 ('({:'+fmt+'}, {:'+fmt+'})').format(cl, cu)])
    
    # attach table
    ax = axs[0]
    ax.set_position([0, 0, f_table, 1])
    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    table = ax.table( 
        cellText = table_content,  
        colLabels = ['$\\bf{Subgroup}$', '$\\bf{N}$', '$\\bf{{{}}}$'.format(label), '$\\bf{CI}$'], 
        loc ='left', cellLoc='center',
        bbox=[0,0, 1, 1])
    # remove gridlines
    for key, cell in table.get_celld().items():
        cell.set_linewidth(0) # remove lines
        cell.fill=False # make transparent for easier post processing
    table.auto_set_column_width(col=list(range(4)))
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    line_x0 = mlines.Line2D([xmin, xmax], [ymin, ymin], color='k')
    line_x1 = mlines.Line2D([xmin, xmax], [ymax*1-1/(n_row+1), ymax*1-1/(n_row+1)], color='k', lw=1)
    ax.add_line(line_x0)
    ax.add_line(line_x1)
    
    # plot errorbar
    ax = axs[1]
    ax.set_position([f_table, 0, 1-f_table, 1-1/(n_row+1)])
    ax.errorbar(X, Y, xerr=[error_lower, error_upper], 
                marker="d", ls='none', color='k',
                capsize=2, elinewidth=.5, markeredgewidth=.5)
    ymin = 0
    ymax = n_row
    ax.set_xlim(xlim)
    ax.set_ylim((ymin, ymax))
    ax.set_xlabel(xlabel)
    # axes lines
    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    line_x0 = mlines.Line2D(xlim, [ymin, ymin], color='k')
    line_x1 = mlines.Line2D(xlim, [ymax, ymax], color='k')
    line_mean = mlines.Line2D([xbar, xbar], [ymin, ymax], color='r', linestyle='dotted')
    ax.add_line(line_x0)
    ax.add_line(line_x1)
    ax.add_line(line_mean)
    return