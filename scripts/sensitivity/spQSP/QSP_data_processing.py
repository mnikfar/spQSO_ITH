# -*- coding: utf-8 -*-
"""
QSP/SPQSP Data pre-processing functions

Created on Wed Jul 29 09:07:38 2020

@author: Chang
"""

import numpy as np

import re
import lxml.etree as ET

import QSP_analysis as qa
import QSP_db as qdb

import matplotlib.pyplot as plt

TAG_QSP_IC_SUCCESS = 'QSP_IC_SUCCESS'
TAG_QSP_IC_FAILURE = 'QSP_IC_FAILURE'
TAG_ABM_SIM_SUCCESS = 'SIMULATION_SUCCESS'

SEC_PER_DAY = 24*3600
MIN_ABM_CC = 1000

#%% data accesing and archiving
def get_lhs_dir(param_dir, group):
    #return param_dir
    return param_dir/'group_{}/'.format(group)

def get_sim_dir(working_dir, group, sample, treat, rep):
    #return working_dir/'sample_{}/rep_{}/'.format(sample,rep)
    #return working_dir/'subject_{}/sample_{}/rep_{}/'.format(sub,sample,rep)
    return working_dir/'group_{}/sample_{}/treatment_{}/rep_{}/'.format(group,sample,treat,rep)
 
# retrieve ids
def get_param_id(db, group, sample):
    condition = 'WHERE group_id = {} AND sample_id = {}'\
        .format(group, sample)
    param_id_list = db.fetch(qdb.TABLE_NAME_PARAM, ['UID'], condition)
    return None if not param_id_list else param_id_list[0][0]

def get_sim_id(db, group, sample, treatment, rep):
    param_id = get_param_id(db, group, sample)
    condition = 'WHERE parameter_id = {} AND treatment_id = {} AND replication_id = {}'\
        .format(param_id, treatment, rep)
    sim_id_list = db.fetch(qdb.TABLE_NAME_SIM, ['UID'], condition)
    return None if not sim_id_list else sim_id_list[0][0]

# save LHS parameter values to database
def archive_param(param_dir, db, group):
    lhs_dir = get_lhs_dir(param_dir, group)
    filename_lhs = [str(f) for f in lhs_dir.glob('param_log.csv')][0]
    header, data = qa.read_csv(filename_lhs)
    cols_param = db.get_colnames(qdb.TABLE_NAME_PARAM)[1:]
    
    for i, row in enumerate(data):
        if get_param_id(db, group, row[0]) is None:
            values = [str(group)] + list(row)
            db.add_entry(qdb.TABLE_NAME_PARAM, cols_param, values)
    db.commit()
    return

# save simulation results to database
def archive_result(working_dir, db, group, sample, treat, rep, windows=['margin','core']):
                
    sim_dir = get_sim_dir(working_dir, group, sample, treat, rep)
    
    # only archive result when initialization is successful
    if not list(sim_dir.glob(TAG_QSP_IC_SUCCESS)):
        print('Failed initialization in:', sim_dir)
        return sim_dir
    
    if not list(sim_dir.glob(TAG_ABM_SIM_SUCCESS)):
        print('Failed abm in:', sim_dir)
        return sim_dir
    
    # Parameter combindation
    param_id = get_param_id(db, group, sample)

    cols_sim = db.get_colnames(qdb.TABLE_NAME_SIM)[1:]
    cols_qsp = db.get_colnames(qdb.TABLE_NAME_QSP)[1:]
    cols_abm = db.get_colnames(qdb.TABLE_NAME_ABM)[1:]
    
    # Simulation    
    sim_id = get_sim_id(db, group, sample, treat, rep)
    if sim_id is None:
        sim_id = db.add_entry(qdb.TABLE_NAME_SIM, cols_sim, [str(param_id), str(treat), str(rep)])
        
        # QSP
        qsp_files = list(sim_dir.glob('QSP*.csv'))
        assert len(qsp_files)==1, 'QSP stats file not found in {}'.format(str(sim_dir))
            
        filename_qsp = str(qsp_files[0])
        
        header, data = qa.read_csv(filename_qsp)
        for row in data:
            values = [str(sim_id)] + list(np.nan_to_num(row.astype(np.float)).astype(str))
            db.add_entry(qdb.TABLE_NAME_QSP, cols_qsp, values)
            
        # ABM core/margin
        for window in windows:
            stats_file_pattern='stats_{}*.csv'.format(window)
            filename_core_abm = [str(f) for f in sim_dir.glob(stats_file_pattern)][0]
            header, data = qa.read_csv(filename_core_abm)
            for row in data:
                values = [str(sim_id), '"{}"'.format(window)] + list(np.nan_to_num(row.astype(np.float)).astype(str))
                db.add_entry(qdb.TABLE_NAME_ABM, cols_abm, values)
                
        db.commit()
    return sim_dir

#%% data pre-processing
######################################################################
# retrieve readouts and preprocess for analysis
######################################################################
# retrieve readouts

READOUT_NAMES = ['d_chg',
                't_to_prog',
                'event_prog',#1 for progression event; 0 for censoring
                'diam_T',
                'tumor_C_end',
                'cent_tcyt',
                'cent_treg',
                'cent_RC_ratio',
                'tumor_tcyt',
                'tumor_treg',
                'tumor_RC_ratio',
                'tumor_C_0',
                'diam_T_0',
                'margin_tcyt',
                'margin_treg',
                'margin_RC_ratio',
                'margin_PDL1',
                'margin_cancer',
                'core_tcyt',
                'core_treg',
                'core_RC_ratio',
                'core_PDL1',
                'core_cancer'
                ]

READOUT_NAMES_H = ['Tumor diameter change',
                'Time to progression',
                'Has progression event',
                'Tumor diameter',
                'Cancer cell number',
                'Tcyt, Blood',
                'Treg, Blood',
                'Treg-cd8 ratio, blood',
                'Tcyt, Tumor',
                'Treg, Tumor',
                'Treg-cd8 ratio, tumor',
                'Cancer cell, Tumor',
                'Tumor diameter',
                'Tcyt, IF',
                'Treg, IF',
                'Treg-cd8 ratio, IF',
                'PDL1, IF',
                'Cancer cell, IF',
                'Tcyt, Core',
                'Treg, Core',
                'Treg-cd8 ratio, Core',
                'PDL1, Core',
                'Cancer cell, Core'
                ]

# return column id from (list of) column names
def get_col(header, namelist):
    if type(namelist) == list:
        idx = [header.index(x) for x in namelist]
        return np.array(idx)
    else:
        idx = header.index(namelist)
        return idx

def abm_scaling(qsp_cc, abm_cc, weight, fraction):
    return qsp_cc/(abm_cc+MIN_ABM_CC) * fraction *(1-weight)/weight

def get_abm_cells(data, header_abm):
    id_abm_cc = get_col(header_abm, ['agentCount.cancerCell.Stem',
                                     'agentCount.cancerCell.Progenitor',
                                     'agentCount.cancerCell.Senescent'])
    id_abm_cd8 =  get_col(header_abm, 
                          ['agentCount.CD8.effector',
                           'agentCount.CD8.cytotoxic',
                           'agentCount.CD8.suppressed'])
    
    id_abm_treg = get_col(header_abm, 'agentCount.Treg.default')
    
    id_PDL1_pos = get_col(header_abm, 'PDL1_pos')
    
    cc = np.sum(data[:, id_abm_cc], axis=1)
    tc = np.sum(data[:, id_abm_cd8], axis=1)
    treg = data[:, id_abm_treg]
    PDL1 = data[:, id_PDL1_pos]

    return cc, tc, treg, PDL1

# calculate tumor diameter change and volume
def get_tum_size(db, sub_id, sample_id, treat_id, rep_id, t_start, weight_qsp):
    sim_id = get_sim_id(db, sub_id, sample_id, treat_id, rep_id)
    # read QSP solution
    condition = 'WHERE sim_id = {} ORDER BY time'.format(sim_id)
    t_qsp = np.array(db.fetch(qdb.TABLE_NAME_QSP, ['time'], condition), dtype=object).flatten()
    row_start = np.where(t_qsp >= t_start)[0][0]
    num_qsp_id_col = 2
    header_qsp = db.get_colnames(qdb.TABLE_NAME_QSP)[num_qsp_id_col:]
    condition = 'WHERE sim_id = {} ORDER BY time'.format(sim_id)
    data_qsp = np.array(db.fetch(qdb.TABLE_NAME_QSP, '*', condition))[:, num_qsp_id_col:].astype(float)
    # column numbers
    id_cc = get_col(header_qsp, 'Tum.C1')
    #id_tcyt_tum = get_col(header_qsp, ['Tum.Teff_1_0','Tum.Teff_exhausted']) 
    id_tcyt_tum = get_col(header_qsp, 'Tum.Teff_1_0')#
    id_tcyt_exh_tum = get_col(header_qsp, 'Tum.Teff_exhausted')#  include exhausted
    id_treg_tum = get_col(header_qsp, 'Tum.Treg')
    id_apc_tum = get_col(header_qsp, ['Tum.APC', 'Tum.mAPC'])
    
    # qsp tumor readouts
    C_qsp = data_qsp[:, id_cc]
    #t_qsp = np.sum(data_qsp[:, id_tcyt_tum], axis=1)
    t_qsp = data_qsp[:, id_tcyt_tum]
    texh_qsp = data_qsp[:, id_tcyt_exh_tum]
    treg_qsp = data_qsp[:, id_treg_tum]
    apc_qsp = np.sum(data_qsp[:, id_apc_tum], axis=1)
        
    cc = C_qsp/weight_qsp
    tc = t_qsp/weight_qsp
    texh = texh_qsp/weight_qsp
    treg = treg_qsp/weight_qsp
    
    d1 = qa.get_tum_diameter(tc+treg+texh, cc, apc_qsp)
    v1 = qa.get_tum_volume(tc+treg+texh, cc, apc_qsp)
    d0 = d1[row_start]
    d_chg = d1/d0 - 1
    return d_chg, v1
    
    
def get_readout(db, sub_id, sample_id, treat_id, rep_id, t_sim_start, weight_qsp, fraction_margin=.1, SLICE_PER_DAY=4):
    
    fraction_core = 1-fraction_margin

    num_readout = len(READOUT_NAMES)
    data = np.zeros(num_readout)
    
    sim_id = get_sim_id(db, sub_id, sample_id, treat_id, rep_id)
    
    # read QSP solution
    condition = 'WHERE sim_id = {} ORDER BY time'.format(sim_id)
    time_qsp = np.array(db.fetch(qdb.TABLE_NAME_QSP, ['time'], condition), dtype=object).flatten()
    row_start = np.where(time_qsp >= t_sim_start)[0][0]
    
    num_qsp_id_col = 2
    header_qsp = db.get_colnames(qdb.TABLE_NAME_QSP)[num_qsp_id_col:]
    condition = 'WHERE sim_id = {} ORDER BY time'.format(sim_id)
    data_qsp = np.array(db.fetch(qdb.TABLE_NAME_QSP, '*', condition))[:, num_qsp_id_col:].astype(float)
    
    # read ABM output
    num_abm_id_col = 3
    header_abm = db.get_colnames(qdb.TABLE_NAME_ABM)[num_abm_id_col:]
    condition = 'WHERE sim_id = {} AND window = "{}" ORDER BY time'.format(sim_id, 'core')
    data_core = np.array(db.fetch(qdb.TABLE_NAME_ABM, '*', condition))[:, num_abm_id_col:].astype(float)

    condition = 'WHERE sim_id = {} AND window = "{}" ORDER BY time'.format(sim_id, 'margin')
    data_margin = np.array(db.fetch(qdb.TABLE_NAME_ABM, '*', condition))[:, num_abm_id_col:].astype(float)
    
    # column numbers
    id_cc = get_col(header_qsp, 'Tum.C1')
    id_tcyt_tum = get_col(header_qsp, ['Tum.Teff_1_0','Tum.Teff_exhausted']) 
    #id_tcyt_tum = get_col(header_qsp, 'Tum.Teff_1_0')# not include exhausted
    id_treg_tum = get_col(header_qsp, 'Tum.Treg')
    id_apc_tum = get_col(header_qsp, ['Tum.APC', 'Tum.mAPC'])
    
    id_tcyt_cent = get_col(header_qsp, 'Cent.Teff_1_0')
    id_treg_cent = get_col(header_qsp, 'Cent.Treg')
    
    # qsp tumor readouts
    C_qsp = data_qsp[:, id_cc]
    t_qsp = np.sum(data_qsp[:, id_tcyt_tum], axis=1)
    #t_qsp = data_qsp[:, id_tcyt_tum]
    treg_qsp = data_qsp[:, id_treg_tum]
    apc_qsp = np.sum(data_qsp[:, id_apc_tum], axis=1)
    all_tumor_qsp = C_qsp+t_qsp+treg_qsp+apc_qsp
    
    # abm readouts
    cc_core, tc_core, treg_core, PDL1_core = get_abm_cells(data_core, header_abm)
    cc_margin, tc_margin, treg_margin, PDL1_margin = get_abm_cells(data_margin, header_abm)
    
    all_core = cc_core + tc_core + treg_core
    all_margin = cc_margin + tc_margin + treg_margin
    
    
    #scaling_core = abm_scaling(C_qsp, cc_core, weight_qsp, fraction_core)
    #scaling_margin = abm_scaling(C_qsp, cc_margin, weight_qsp, fraction_margin)

    #cc = C_qsp + scaling_core*cc_core + scaling_margin*cc_margin
    cc = C_qsp/weight_qsp
    #tc = t_qsp + scaling_core*tc_core + scaling_margin*tc_margin
    tc = t_qsp/weight_qsp
    #treg = treg_qsp + scaling_core*treg_core + scaling_margin*treg_margin
    treg = treg_qsp/weight_qsp
    
    d1 = qa.get_tum_diameter(tc+treg, cc, apc_qsp)
    v1 = qa.get_tum_volume(tc+treg, cc, apc_qsp)
    d0 = d1[row_start]
    d_chg = d1/d0 - 1


    time_max_day = max(time_qsp) + 1
    data[get_col(READOUT_NAMES,'d_chg')] = d_chg[-1]
    #argmax(d_chg>.2): index of first time point where diameter change is larger than 20%
    prog_event_true = sum(d_chg>.2) > 0
    data[get_col(READOUT_NAMES,'t_to_prog')] = ((time_qsp[np.argmax(d_chg>.2)] if prog_event_true else time_max_day)-t_sim_start)/SLICE_PER_DAY
    data[get_col(READOUT_NAMES,'event_prog')] = prog_event_true
    data[get_col(READOUT_NAMES,'diam_T')] = d1[-1]
    data[get_col(READOUT_NAMES,'tumor_C_end')] = cc[-1]
    data[get_col(READOUT_NAMES,'cent_tcyt')] = data_qsp[row_start, id_tcyt_cent]
    data[get_col(READOUT_NAMES,'tumor_tcyt')] = t_qsp[row_start]
    data[get_col(READOUT_NAMES,'cent_treg')] = data_qsp[row_start, id_treg_cent]
    data[get_col(READOUT_NAMES,'tumor_treg')] = treg_qsp[row_start]
    data[get_col(READOUT_NAMES,'cent_RC_ratio')] = data_qsp[row_start, id_treg_cent]/data_qsp[row_start, id_tcyt_cent]
    data[get_col(READOUT_NAMES,'tumor_RC_ratio')] = treg_qsp[row_start]/t_qsp[row_start]
    data[get_col(READOUT_NAMES,'tumor_C_0')] = cc[row_start]
    data[get_col(READOUT_NAMES,'diam_T_0')] = d1[row_start]
    
    data[get_col(READOUT_NAMES,'margin_tcyt')] = tc_margin[row_start]
    data[get_col(READOUT_NAMES,'margin_treg')] = treg_margin[row_start]
    data[get_col(READOUT_NAMES,'margin_cancer')] = cc_margin[row_start]
    data[get_col(READOUT_NAMES,'margin_RC_ratio')] = treg_margin[row_start]/tc_margin[row_start]
    data[get_col(READOUT_NAMES,'margin_PDL1')] = PDL1_margin[row_start]/all_margin[row_start]
    
    data[get_col(READOUT_NAMES,'core_tcyt')] = tc_core[row_start]
    data[get_col(READOUT_NAMES,'core_treg')] =  treg_core[row_start]
    data[get_col(READOUT_NAMES,'core_cancer')] = cc_core[row_start]
    data[get_col(READOUT_NAMES,'core_RC_ratio')] = treg_core[row_start]/tc_core[row_start]
    data[get_col(READOUT_NAMES,'core_PDL1')] = PDL1_core[row_start]/all_core[row_start]
    
    if tc_margin[row_start] == 0 or tc_core[row_start]==0:
        raise ValueError('0 tcyt in ratio')
        
    return data, d_chg, d1, v1

# get parameter values as a dictionary
def get_param_dict(param_file):
    
    master_value_pattern = re.compile('^\s*\{.*\}\s*$')
    value_range_pattern = re.compile('^\s*\[.*,.*\]\s*$')
    ATTRIB_STAGE_NAME = 'stage'
    STAGE_TAG_PRE = 'pre'
    STAGE_TAG_POST = 'post'

    param_dict = {}
    param_dict_const = {}
    param_dict_group = {}
    param_dict_treat = {}
    param_dict['const'] = param_dict_const
    param_dict['group'] = param_dict_group
    param_dict['treat'] = param_dict_treat

    tree = ET.parse(param_file)
    root = tree.getroot()

    for elem in root.iter():
        if not len(elem):
            path = tree.getelementpath(elem)
            #param_key = ('.').join(path.split('/'))
            param_key = path
            if value_range_pattern.match(elem.text):
                pass
            elif master_value_pattern.match(elem.text):
                stage = elem.attrib[ATTRIB_STAGE_NAME]
                s = elem.text.split('}')[0].split('{')[1]
                a = np.asarray(s.split(',')).astype('float')
                if stage == STAGE_TAG_PRE:
                    param_dict_group[param_key] = a
                elif stage == STAGE_TAG_POST:
                    param_dict_treat[param_key] = a
                else:
                    raise ValueError('Unknown master stage type: {}'.format(stage))
            else:
                try:
                    param_val = float(elem.text)
                    param_dict_const[param_key] = param_val
                except:
                    pass
    return param_dict

def get_param_val(db, name, param_dict, group, sample, treat):
    header_lhs, lhs = get_LHS(db, group)
    val = 0
    if name in header_lhs:
        idx_param = header_lhs.index(name)
        val = lhs[sample-1, idx_param]
    else:
        if name in param_dict['const']:
            val =  param_dict['const'][name]
        elif name in param_dict['group']:
            val = param_dict['group'][name][group-1]
        elif name in param_dict['treat']:
            val =  param_dict['treat'][name][treat-1]
        else:
            raise ValueError('Unknown parameter: {}'.format(name))
    return val

# get LHS parameter name and values
def get_LHS(db, group):
    num_lhs_col = 2
    header_lhs = db.get_colnames(qdb.TABLE_NAME_PARAM)[num_lhs_col:]
    condition = 'WHERE group_id = {} ORDER BY sample_id'.format(group)
    lhs = np.array(db.fetch(qdb.TABLE_NAME_PARAM, '*', condition))[:, num_lhs_col:].astype(float)
    return header_lhs, lhs

# get single output time series
def get_qsp_time_series(db, sub_id, sample_id, treat_id, rep_id, colname):
    sim_id = get_sim_id(db, sub_id, sample_id, treat_id, rep_id)
    data_qsp = get_qsp_time_series_by_sim_id(db, sim_id, [colname])
    return data_qsp.flatten()

# get single output time series
# colname can be a substring of column name. 
# In this case, the sum of all matching cols is reported  
def get_abm_time_series(db, sub_id, sample_id, treat_id, rep_id, window, colname):
    sim_id = get_sim_id(db, sub_id, sample_id, treat_id, rep_id)
    abm_cols = db.get_colnames(qdb.TABLE_NAME_ABM)
    col_list = [s for s in abm_cols if colname in s]
    #print(col_list)
    data_abm = get_abm_time_series_by_sim_id(db, sim_id, window, col_list)
    #print(data_abm.shape)
    return np.sum(data_abm, 1)

# increment between two time points
def get_abm_time_series_inc(db, sub_id, sample_id, treat_id, rep_id, window, colname):
    data_abm_flatten = get_abm_time_series(db, sub_id, sample_id, treat_id, rep_id, window, colname)
    res = np.zeros(data_abm_flatten.shape)
    res[1:] = data_abm_flatten[1:] - data_abm_flatten[:-1]
    return res

# column name need to be in a list
def get_qsp_time_series_by_sim_id(db, sim_id, cols):
    condition = 'WHERE sim_id = {} ORDER BY time'.format(sim_id)
    data_qsp = np.array(db.fetch(qdb.TABLE_NAME_QSP, 
                                 cols, condition)).astype(float)
    return data_qsp

# column name need to be in a list
def get_abm_time_series_by_sim_id(db, sim_id, window, cols):
    condition = 'WHERE sim_id = {} and window = "{}" ORDER BY time'.format(sim_id, window)
    #print(condition)
    data_abm = np.array(db.fetch(qdb.TABLE_NAME_ABM, 
                                 cols, condition)).astype(float)
    return data_abm
    
def get_result(data, colstr):
    return data[:,READOUT_NAMES.index(colstr)]

#%% plotting
    
def plot_general(db, params, group, treat, samples, reps, colors, alpha):
    
    NIVO_START_DAY = get_param_val(db, 'QSP/init_value/Parameter/t_init_Nivo', params, 1, 1, 1)
    SLICE_PER_DAY = SEC_PER_DAY / get_param_val(db, 'ABM/Environment/SecPerSlice', params, 1, 1, 1)
    t_nivo_start = NIVO_START_DAY*SLICE_PER_DAY

    use_auto_ylim = False
    V2A = 0.005+0.01-2*0.0025
    V_C = get_param_val(db, 'QSP/init_value/Parameter/vol_cent', params, 1, 1, 1) * 1e6 #mm^3

    sim_id_0 = db.fetch(qdb.TABLE_NAME_ABM, ['UID'], '')[0][0]
    t_qsp_day = get_qsp_time_series_by_sim_id(db, sim_id_0, ['time']).flatten()/SLICE_PER_DAY

    # direct results
    var_list = [# row 1
                {'loc': (0, 0), 'scale': 1, 'ylim': (1e1, 1e11), 'col':'Tum.C1', 'title': 'Cancer cell'},
                {'loc': (0, 3), 'scale': 1e12, 'ylim': (1e2, 1e8), 'col': 'Tum.D1_0', 'title': 'Tumor antigen, pM'},
                {'loc': (0, 4), 'scale': 1, 'ylim': (1e1, 1e10), 'col': 'Tum.mAPC', 'title': 'Tumor mAPC' },
                # row 2
                {'loc': (1, 0), 'scale': 1, 'ylim': (1e1, 1e10), 'col': 'Tum.Teff_1_0', 'title': 'Tumor Teff' },
                {'loc': (1, 1), 'scale': 1, 'ylim': (1e1, 1e10), 'col': 'Tum.Treg', 'title': 'Tumor Treg' },
                {'loc': (1, 2), 'scale': 1, 'ylim': (1e1, 1e10), 'col': 'Tum.Teff_exhausted', 'title': 'Tumor exhuasted Teff' },
                # row 4
                {'loc': (3, 0), 'scale': 1, 'ylim': (1e0, 1e10), 'col': 'LN.aT_1_0', 'title': 'LN activated Teff' },
                {'loc': (3, 1), 'scale': 1, 'ylim': (1e4, 1e14), 'col': 'LN.Teff_1_0', 'title': 'LN mature Teff'},
                {'loc': (3, 2), 'scale': 1, 'ylim': (1e2, 1e14), 'col': 'Cent.Teff_1_0', 'title': 'Blood Teff' },
                {'loc': (3, 3), 'scale': 1/V_C, 'ylim': (1e-3, 1e5), 'col': 'Cent.Teff_1_0', 'title': 'Blood Teff, $mm^{-3}$'},
                {'loc': (3, 4), 'scale': 1, 'ylim': (1e6, 1e16), 'col': 'Peri.Teff_1_0', 'title':  'Peripheral  Teff'},
                # row 5
                {'loc': (4, 0), 'scale': 1, 'ylim': (1e1, 1e10), 'col': 'LN.aTreg_CD4', 'title': 'LN activated Treg' },
                {'loc': (4, 1), 'scale': 1, 'ylim': (1e2, 1e12), 'col': 'LN.Treg', 'title': 'LN mature Treg'},
                {'loc': (4, 2), 'scale': 1, 'ylim': (1e0, 1e10), 'col': 'Cent.Treg', 'title': 'Blood Treg' },
                {'loc': (4, 3), 'scale': 1/V_C, 'ylim': (1e-3, 1e3), 'col': 'Cent.Treg', 'title': 'Blood Treg, $mm^{-3}$'},
                {'loc': (4, 4), 'scale': 1, 'ylim': (1e6, 1e16), 'col': 'Peri.Treg', 'title':  'Peripheral Treg'}
                ]
    
    # figure
    #fig, axs = plt.subplots(5,5, constrained_layout=True)
    fig, axs = plt.subplots(5,5)
    
    for s in samples:
        sim_id = get_sim_id(db, group, s, treat, 1) # None if no such simulation archived
        if sim_id is not None:
            for i, var in enumerate(var_list):
                ax = axs[var['loc']]
                y = get_qsp_time_series(db, group, s, treat, 1, var['col'])
                ax.plot(t_qsp_day, y*var['scale'],  alpha=alpha, color = colors[s-1])
                ax.set_title(var['title'])
                ax.set_xlabel('t, day')
                ax.set_yscale('log')
                if var['ylim'] and not use_auto_ylim:
                    ax.set_ylim(*var['ylim'])
            
            Ccell = get_qsp_time_series(db, group, s, treat, 1, 'Tum.C1')
            Ccell[Ccell<1]=1
            
            weight_qsp = get_param_val(db, 'QSP/simulation/weight_qsp',  params, group, s, treat)
            
            T1_tum = get_qsp_time_series(db, group, s, treat, 1, 'Tum.Teff_1_0')
            Tr_tum = get_qsp_time_series(db, group, s, treat, 1, 'Tum.Treg')
            Te_tum = get_qsp_time_series(db, group, s, treat, 1, 'Tum.Teff_exhausted')
            d_chg_time, vol_time = get_tum_size(db, group, s, treat, 1, t_nivo_start, weight_qsp)
            
            # Row 1
            # spider
            ax = axs[0, 1]
            line = ax.plot(t_qsp_day-NIVO_START_DAY, d_chg_time, alpha=alpha, color = colors[s-1])
            ax.set_xlim(0, t_qsp_day[-1]-NIVO_START_DAY)
            ax.set_ylim(-1, 1.5)
            ax.set_title('Tumor diameter change')
            
            # vol
            ax = axs[0, 2]
            line = ax.plot(t_qsp_day, vol_time, alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_title('Tumor volume, $mm^{3}$')
            
            # Row 2
            #Teff/all
            ax = axs[1, 3]
            ax.plot(t_qsp_day, T1_tum/(T1_tum+Te_tum+Tr_tum+Ccell)*100, alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(1e-3, 100)
            ax.set_title('Tumor Teff / Total cell, %')
            
            #Teff/(Texh+Teff)
            ax = axs[1, 4]
            ax.plot(t_qsp_day, T1_tum/(Te_tum+T1_tum)*100, alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(1e-2, 200)
            ax.set_title('Tumor Teff / Teff_all, %')
            
            # Row 3
            #Teff n/V
            ax = axs[2, 0]
            Teff_dv = (T1_tum+Te_tum)/vol_time
            ax.plot(t_qsp_day, Teff_dv, alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(.1, 1e8)
            ax.set_title('Tumor Teff volume density, $mm^{-3}$')
            
            #Treg n/V
            ax = axs[2, 1]
            Treg_dv = (Tr_tum)/vol_time
            ax.plot(t_qsp_day, Treg_dv, alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(.1, 1e8)
            ax.set_title('Tumor Treg  volume density, $mm^{-3}$')
            
            #Teff  n/A
            ax = axs[2, 2]
            ax.plot(t_qsp_day, Teff_dv*V2A, alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(.1, 1e5)
            ax.set_title('Tumor Teff area density, $mm^{-2}$')
            
            #Treg n/A
            ax = axs[2, 3]
            ax.plot(t_qsp_day, Treg_dv*V2A, alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(.1, 1e5)
            ax.set_title('Tumor Treg area density, $mm^{-2}$')
            
            #Treg/Teff
            ax = axs[2, 4]
            ax.plot(t_qsp_day, 100 * Tr_tum/(T1_tum+Te_tum), alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(1e-3, 1e3)
            ax.set_title('Tumor Treg/Teff, %')
            
            # Row 4
            # Row 5
            
    num_row = 5
    num_col = 5
    for r in range(num_col):
        axs[num_row-1,0].set_xlabel('t, day')
    
    fig.set_size_inches(20,15)
    fig.tight_layout()
    
    return fig

##%% plot priming
def plot_priming(db, params, group, treat, samples, reps, colors, alpha):
    
    SLICE_PER_DAY = SEC_PER_DAY / get_param_val(db, 'ABM/Environment/SecPerSlice', params, 1, 1, 1)
    AVO = get_param_val(db, 'QSP/init_value/Parameter/AvogadroN', params, 1, 1, 1)

    sim_id_0 = db.fetch(qdb.TABLE_NAME_ABM, ['UID'], '')[0][0]
    t_qsp_day = get_qsp_time_series_by_sim_id(db, sim_id_0, ['time']).flatten()/SLICE_PER_DAY
    
    
    use_auto_ylim = False
    var_list = [# row 1
                {'loc': (0, 0), 'scale': 1, 'ylim': (), 'col':'Tum.Ckine_Mat', 'title': 'Tum.Ckine_Mat'},
                {'loc': (0, 1), 'scale': 1e12, 'ylim': (1e0, 1e8), 'col': 'Tum.D1_0', 'title': 'Tumor antigen, pM'},
                {'loc': (0, 2), 'scale': 1, 'ylim': (), 'col': 'Tum.Ckine_Mat', 'title': 'Tumor Ckine, nM'},
                {'loc': (0, 3), 'scale': 1, 'ylim': (), 'col': 'Tum.mAPC', 'title': 'Tum.mAPC'},
                {'loc': (0, 4), 'scale': 1, 'ylim': (), 'col': 'LN.mAPC', 'title': 'LN.mAPC'},
                {'loc': (1, 0), 'scale': 1, 'ylim': (), 'col': 'LN.IL2', 'title': 'LN.IL2, nM'},
                {'loc': (1, 3), 'scale': 1, 'ylim': (), 'col': 'Tum.APC', 'title': 'Tum.APC'},
                {'loc': (1, 4), 'scale': AVO, 'ylim': (1e2, 1e8), 'col': 'mAPC_Surf.p1_0_M1', 'title': 'mAPC_Surf.p1_0_M'},
                {'loc': (2, 4), 'scale': 1, 'ylim': (), 'col': 'LN.nT_CD8', 'title': 'LN.nT_CD8'},
                {'loc': (3, 1), 'scale': 1, 'ylim': (), 'col': 'LN.aT_1_0', 'title': 'LN.aTeff'},
                {'loc': (3, 2), 'scale': 1, 'ylim': (), 'col': 'LN.Teff_1_0', 'title': 'LN.Teff'},
                {'loc': (3, 3), 'scale': 1, 'ylim': (), 'col': 'Cent.Teff_1_0', 'title': 'Cent.Teff'},
                {'loc': (3, 4), 'scale': 1, 'ylim': (), 'col': 'Peri.Teff_1_0', 'title': 'Peri.Teff'},
                {'loc': (4, 1), 'scale': 1, 'ylim': (), 'col': 'LN.aTreg_CD4', 'title': 'LN.aTreg'},
                {'loc': (4, 2), 'scale': 1, 'ylim': (), 'col': 'LN.Treg', 'title': 'LN.Treg'},
                {'loc': (4, 3), 'scale': 1, 'ylim': (), 'col': 'Cent.Treg', 'title': 'Cent.Treg'},
                {'loc': (4, 4), 'scale': 1, 'ylim': (), 'col': 'Peri.Treg', 'title': 'Peri.Treg'}
                ]
    fig, axs = plt.subplots(5,5)
    
    for s in samples:
        sim_id = get_sim_id(db, group, s, treat, 1) # None if no such simulation archived
        if sim_id is not None:
            # parameter values
            n_sites_mAPC = get_param_val(db, 'QSP/init_value/Parameter/n_sites_mAPC',  params, group, s, treat) # dimensionless
            nT_CD8_diver = get_param_val(db, 'QSP/init_value/Parameter/nT_CD8_diver',  params, group, s, treat) # dimensionless
            n_clone_p1_0 = get_param_val(db, 'QSP/init_value/Parameter/n_clone_p1_0',  params, group, s, treat) # dimensionless
            K_p_M = get_param_val(db, 'QSP/init_value/Parameter/K_p_M',  params, group, s, treat) # molecule
            k_nTCD8_mAPC = get_param_val(db, 'QSP/init_value/Parameter/k_nTCD8_mAPC',  params, group, s, treat) # day^-1
            n_sites_APC = get_param_val(db, 'QSP/init_value/Parameter/n_sites_APC',  params, group, s, treat) # dimensionless
            nT_CD4_diver = get_param_val(db, 'QSP/init_value/Parameter/nT_CD4_diver',  params, group, s, treat) # dimensionless
            n_clone_Treg = get_param_val(db, 'QSP/init_value/Parameter/n_clone_Treg',  params, group, s, treat) # dimensionless
            k_nTCD4_APC = get_param_val(db, 'QSP/init_value/Parameter/k_nTCD4_APC',  params, group, s, treat) # day^-1
            #k_aT_prolif =  get_param_val(db, 'QSP/init_value/Parameter/k_aT_prolif',  params, group, s, treat) # day^-1
            #n_aT_prolif_TCR = get_param_val(db, 'QSP/init_value/Parameter/n_aT_prolif_TCR',  params, group, s, treat) # dimensionless
            #n_aT_prolif_Coestim =get_param_val(db, 'QSP/init_value/Parameter/n_aT_prolif_Coestim',  params, group, s, treat) # dimensionless
            n_aT_prolif_IL20 = get_param_val(db, 'QSP/init_value/Parameter/n_aT_prolif_IL20',  params, group, s, treat) # dimensionless
            n_Treg_prolif_IL20 = get_param_val(db, 'QSP/init_value/Parameter/n_Treg_prolif_IL20',  params, group, s, treat) # dimensionless
            f_vol_LN =get_param_val(db, 'QSP/init_value/Parameter/f_vol_LN',  params, group, s, treat) # dimensionless
            K_IL2_Teff = get_param_val(db, 'QSP/init_value/Parameter/K_IL2_Teff',  params, group, s, treat) # nmole/L
            K_IL2_Treg = get_param_val(db, 'QSP/init_value/Parameter/K_IL2_Treg',  params, group, s, treat) # nmole/L
    
            for i, var in enumerate(var_list):
                ax = axs[var['loc']]
                y = get_qsp_time_series(db, group, s, treat, 1, var['col'])
                ax.plot(t_qsp_day, y*var['scale'],  alpha=alpha, color = colors[s-1])
                ax.set_title(var['title'])
                ax.set_xlabel('t, day')
                ax.set_yscale('log')
                if var['ylim'] and not use_auto_ylim:
                    ax.set_ylim(*var['ylim'])
            
            # n_prolif_IL2
            ax = axs[1, 1]
            IL2 = get_qsp_time_series(db, group, s, treat, 1, 'LN.IL2')
            n_aT_prolif_IL2 = n_aT_prolif_IL20 * (IL2 / f_vol_LN) / (IL2 / f_vol_LN + K_IL2_Teff) + 0.001
            #n_aT_prolif = n_aT_prolif_IL2 + n_aT_prolif_TCR + n_aT_prolif_Coestim
            ax.plot(t_qsp_day, n_aT_prolif_IL2,  alpha=alpha, color = colors[s-1])
            ax.set_ylim(0, n_aT_prolif_IL20+1)
            ax.set_title('n_aT_prolif_IL2')
            ax.set_xlabel('t, day')
            
            # n_Treg_prolif_IL2
            ax = axs[1, 2]
            n_Treg_prolif_IL2 = n_Treg_prolif_IL20 * (IL2 / f_vol_LN) / (IL2 / f_vol_LN + K_IL2_Treg)
            #n_aT_prolif = n_aT_prolif_IL2 + n_aT_prolif_TCR + n_aT_prolif_Coestim
            ax.plot(t_qsp_day, n_Treg_prolif_IL2,  alpha=alpha, color = colors[s-1])
            ax.set_ylim(0, n_Treg_prolif_IL20+1)
            ax.set_title('n_Treg_prolif_IL2')
            ax.set_xlabel('t, day')
            
            # H_P1
            ax = axs[2, 0]
            p1MHC = get_qsp_time_series(db, group, s, treat, 1, 'mAPC_Surf.p1_0_M')
            H_P1 = p1MHC / n_clone_p1_0 / (p1MHC / n_clone_p1_0 + K_p_M)
            ax.plot(t_qsp_day, H_P1,  alpha=alpha, color = colors[s-1])
            ax.set_title('H_P1')
            ax.set_ylim(0, 1.05)
            ax.set_xlabel('t, day')
                
            # H_mAPC
            ax = axs[2, 1]
            mAPC = get_qsp_time_series(db,group, s, treat, 1, 'LN.mAPC')
            nT_CD8 = get_qsp_time_series(db, group, s, treat, 1, 'LN.nT_CD8')
            H_mAPC = n_sites_mAPC * mAPC / (n_sites_mAPC * mAPC + nT_CD8 / nT_CD8_diver)
            ax.plot(t_qsp_day, H_mAPC,  alpha=alpha, color = colors[s-1])
            ax.set_title('H_mAPC')
            ax.set_ylim(0, 1.05)
            ax.set_xlabel('t, day')
            
            # H_P0
            ax = axs[2, 2]
            p0MHC = get_qsp_time_series(db, group, s, treat, 1, 'mAPC_Surf.cpt_M')
            H_P0 = p0MHC / n_clone_Treg / (p0MHC / n_clone_Treg + K_p_M)
            ax.plot(t_qsp_day, H_P0,  alpha=alpha, color = colors[s-1])
            ax.set_title('H_P0')
            ax.set_ylim(0, 1.05)
            ax.set_xlabel('t, day')
                
            # H_APC
            ax = axs[2, 3]
            APC = get_qsp_time_series(db,group, s, treat, 1, 'LN.APC')
            nT_CD4 = get_qsp_time_series(db, group, s, treat, 1, 'LN.nT_CD4')
            H_APC = n_sites_APC * APC / (n_sites_APC * APC + nT_CD4 / nT_CD4_diver)
            ax.plot(t_qsp_day, H_APC,  alpha=alpha, color = colors[s-1])
            ax.set_title('H_APC')
            ax.set_ylim(0, 1.05)
            ax.set_xlabel('t, day')
            
           
            # nCD8 activation
            ax = axs[3, 0]
            nTeff = get_qsp_time_series(db, group, s, treat, 1, 'LN.nT_CD8')
            act_Teff = k_nTCD8_mAPC* n_clone_p1_0*(nTeff/nT_CD8_diver)*H_mAPC*H_P1
            ax.plot(t_qsp_day, act_Teff,  alpha=alpha, color = colors[s-1])
            ax.set_title('nCD8 activation')
            ax.set_xlabel('t, day')
            
            # nCD4 activation
            ax = axs[4, 0]
            nTreg = get_qsp_time_series(db, group, s, treat, 1, 'LN.nT_CD4')
            act_Treg = k_nTCD4_APC* n_clone_Treg*(nTreg/nT_CD4_diver)*H_APC*H_P0
            ax.plot(t_qsp_day, act_Treg,  alpha=alpha, color = colors[s-1])
            ax.set_title('nCD4 activation')
            ax.set_xlabel('t, day')
    
    
    fig.set_size_inches(20,15)
    fig.tight_layout()
    return fig

##%% killing mechanism
def plot_killing(db, params, group, treat, samples, reps, colors, alpha):    
    fig, axs = plt.subplots(4,4)
    
    SLICE_PER_DAY = SEC_PER_DAY / get_param_val(db, 'ABM/Environment/SecPerSlice', params, 1, 1, 1)
    sim_id_0 = db.fetch(qdb.TABLE_NAME_ABM, ['UID'], '')[0][0]
    t_qsp_day = get_qsp_time_series_by_sim_id(db, sim_id_0, ['time']).flatten()/SLICE_PER_DAY
    
    use_auto_ylim = False
    
    var_list = [{'loc': (1, 0), 'scale': 1, 'ylim': (1e0, 1e6), 'col':'Tum.C1_PDL1_Teff_PD1', 'title': 'Tum.C1_PDL1_Teff_PD1'},
                {'loc': (1, 1), 'scale': 1, 'ylim': (1e0, 1e6), 'col':'Tum.C1_PDL2_Teff_PD1', 'title': 'Tum.C1_PDL2_Teff_PD1'},
                {'loc': (1, 2), 'scale': 1, 'ylim': (1e0, 1e6), 'col':'Tum.C1_PDL1_syn', 'title':'Tum.C1_PDL1_syn' },
                {'loc': (1, 3), 'scale': 1, 'ylim': (1e0, 1e6), 'col':'Tum.C1_PDL2_syn', 'title': 'Tum.C1_PDL2_syn'},
                {'loc': (2, 0), 'scale': 1, 'ylim': (1e0, 1e6), 'col':'Tum.Teff_PD1_syn', 'title': 'Tum.Teff_PD1_syn'},
                {'loc': (2, 1), 'scale': 1, 'ylim': (1e0, 1e8), 'col':'Tum.Teff_PD1_Nivo_syn', 'title':'Tum.Teff_PD1_Nivo_syn', },
                {'loc': (2, 2), 'scale': 1, 'ylim': (1e0, 1e5), 'col':'Tum.Teff_PD1_Nivo_PD1_syn', 'title': 'Tum.Teff_PD1_Nivo_PD1_syn'},
                {'loc': (2, 3), 'scale': 1e9, 'ylim': (1e1, 1e3), 'col':'Tum.Nivo', 'title': 'Tum.Nivo, nM'},
                {'loc': (3, 0), 'scale': 1, 'ylim': (1e0, 1e6), 'col': 'Tum.Teff_PD1', 'title':'Tum.Teff_PD1' },
                {'loc': (3, 1), 'scale': 1, 'ylim': (1e0, 1e6), 'col':'Tum.Teff_PD1_Nivo', 'title': 'Tum.Teff_PD1_Nivo'}]
    
    
    for s in samples:
        sim_id = get_sim_id(db, group, s, treat, 1) # None if no such simulation archived
        if sim_id is not None:
            
            k_C_death_by_T = get_param_val(db, 'QSP/init_value/Parameter/k_C_death_by_T',  params, group, s, treat) # 1/day
            k_Teff_death_by_C = get_param_val(db, 'QSP/init_value/Parameter/k_Teff_death_by_C',  params, group, s, treat)# 1/day
            k_Teff_inhibBy_Treg = get_param_val(db, 'QSP/init_value/Parameter/k_Teff_inhibBy_Treg',  params, group, s, treat) #1/day
            K_C1_PDLX_Teff_PD1 = get_param_val(db, 'QSP/init_value/Parameter/K_C1_PDLX_Teff_PD1',  params, group, s, treat) # item
            n_PD1_PDLX = get_param_val(db, 'QSP/init_value/Parameter/n_PD1_PDLX',  params, group, s, treat) # dimensionless
            k_C_death = get_param_val(db, 'QSP/init_value/Parameter/k_C_death',  params, group, s, treat) # 1/day
            k_C_growth = get_param_val(db, 'QSP/init_value/Parameter/k_C_growth',  params, group, s, treat) # 1/day
            vol_tum_max = get_param_val(db, 'QSP/init_value/Parameter/vol_tum_max',  params, group, s, treat) # L
            D_C = get_param_val(db, 'QSP/init_value/Parameter/D_C',  params, group, s, treat) # micron
            f_vol_tum = get_param_val(db, 'QSP/init_value/Parameter/f_vol_tum',  params, group, s, treat) # L
            
            
            for i, var in enumerate(var_list):
                ax = axs[var['loc']]
                y = get_qsp_time_series(db, group, s, treat, 1, var['col'])
                ax.plot(t_qsp_day, y*var['scale'],  alpha=alpha, color = colors[s-1])
                ax.set_title(var['title'])
                ax.set_xlabel('t, day')
                ax.set_yscale('log')
                if var['ylim'] and not use_auto_ylim:
                    ax.set_ylim(*var['ylim'])
            
            vol_C = 3.1416/6*D_C**3 * 1e-15 # L
            K_C_MAX = vol_tum_max/vol_C*(1-f_vol_tum)
            tum_C = get_qsp_time_series(db, group, s, treat, 1,'Tum.C1') # cell
            tum_Treg = get_qsp_time_series(db, group, s, treat, 1,'Tum.Treg') # cell
            tum_Teff = get_qsp_time_series(db, group, s, treat, 1,'Tum.Teff_1_0') # cell
            tum_T_exh = get_qsp_time_series(db, group, s, treat, 1,'Tum.Teff_exhausted') #cell
            tum_all = tum_C + tum_Teff + tum_T_exh + tum_Treg # exhausted T cell accounted
            PDLX_PD1 = get_qsp_time_series(db, group, s, treat, 1,'Tum.C1_PDL1_Teff_PD1') + \
                     get_qsp_time_series(db, group, s, treat, 1, 'Tum.C1_PDL2_Teff_PD1') # molecule
            
            
            H_PDLX_PD1_n = f_hill(PDLX_PD1, K_C1_PDLX_Teff_PD1, n_PD1_PDLX)
            H_PDLX_PD1_1 = f_hill(PDLX_PD1, K_C1_PDLX_Teff_PD1, 1)
            
            dcdt_g = k_C_growth * (1-tum_C/K_C_MAX)
            dcdt_d = k_C_death
            dcdt_k = k_C_death_by_T * tum_Teff/tum_all *(1-H_PDLX_PD1_n)
            
            #dC/Cdt: all
            ax = axs[0, 0]
            dcdt = dcdt_g - dcdt_d - dcdt_k
            ax.plot(t_qsp_day, dcdt,  alpha=alpha, color = colors[s-1])
            ax.set_yscale('symlog', linthreshy=1e-3)
            ax.set_ylim(-1e1, 1e-1)
            #ax.set_ylim(1e-3, 1e1)
            ax.set_title('dC/Cdt, 1/day')
            ax.set_xlabel('t, day')
            
            #dC/Cdt: killing
            ax = axs[0, 1]
            
            ax.plot(t_qsp_day, dcdt_k,  alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(1e-4, 1e1)
            ax.set_title('dC/Cdt, killed, 1/day')
            ax.set_xlabel('t, day')
            
            # H 
            ax = axs[0, 2]
            ax.plot(t_qsp_day, H_PDLX_PD1_n,  alpha=alpha, color = colors[s-1])
            ax.set_yscale('linear')
            ax.set_ylim(0.1, 1.1)
            ax.set_title('H_PDL1_PD1')
            ax.set_xlabel('t, day')
            
            ax = axs[0, 3]
            ax.plot(t_qsp_day, tum_Teff/tum_all,  alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(1e-5, 1.1e0)
            ax.set_title('Teff/Total cell')
            ax.set_xlabel('t, day')
            
            exh_by_C = k_Teff_death_by_C * tum_C/tum_all * H_PDLX_PD1_1
            ax = axs[3, 2]
            ax.plot(t_qsp_day, exh_by_C,  alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(1e-4, 1e0)
            ax.set_title('Exhaustion by C, 1/day')
            ax.set_xlabel('t, day')
            
            inhib_by_reg = k_Teff_inhibBy_Treg *tum_Treg/tum_all*(1+H_PDLX_PD1_n)
            ax = axs[3, 3]
            ax.plot(t_qsp_day, inhib_by_reg,  alpha=alpha, color = colors[s-1])
            ax.set_yscale('log')
            ax.set_ylim(1e-4, 1e4)
            ax.set_title('Exhaustion by Treg, 1/day')
            ax.set_xlabel('t, day')
            
    fig.set_size_inches(20,15)
    fig.tight_layout()
    return fig

def plot_abm(db, params, group, treat, samples, reps, colors, alpha, windows=['margin', 'core']):

    SLICE_PER_DAY = SEC_PER_DAY / get_param_val(db, 'ABM/Environment/SecPerSlice', params, 1, 1, 1)
    sim_id_0 = db.fetch(qdb.TABLE_NAME_ABM, ['UID'], '')[0][0]
    t_qsp_day = get_qsp_time_series_by_sim_id(db, sim_id_0, ['time']).flatten()/SLICE_PER_DAY
    DPS = t_qsp_day[1]-t_qsp_day[0]
    
    s = samples[0]
    size_X = get_param_val(db, 'ABM/Environment/Tumor/XSize',  params, group, s, treat)
    size_Y = get_param_val(db, 'ABM/Environment/Tumor/XSize',  params, group, s, treat) 
    size_Z = get_param_val(db, 'ABM/Environment/Tumor/XSize',  params, group, s, treat) 
    len_voxel = get_param_val(db, 'ABM/Environment/Tumor/VoxelSize',  params, group, s, treat)
    ABM_Vol = size_X * size_Y * size_Z * len_voxel**3 * 1e-9
    
    #print('DPS: ', DPS)
    #AVO = get_param_val(db, 'QSP/init_value/Parameter/AvogadroN', params, 1, 1, 1)
    use_auto_ylim = False
    num_window = len(windows)
    # direct results
    var_list = [# row 1
                {'loc': (0, 0), 'scale': 1/ABM_Vol, 'ylim': (), 'type': 'raw','col':'agentCount.cancerCell', 'title': 'Cancer cell, '+r'$1/mm^{3}$', },
                {'loc': (0, 1), 'scale': 1/ABM_Vol, 'ylim': (), 'type': 'raw','col':'agentCount.cancerCell.Stem', 'title': 'Cancer cell, Stem, '+r'$1/mm^{3}$', },
                {'loc': (0, 2), 'scale': 1/ABM_Vol, 'ylim': (), 'type': 'raw', 'col':'agentCount.cancerCell.Progenitor', 'title': 'Cancer cell, Progenitor, '+r'$1/mm^{3}$'},
                {'loc': (0, 3), 'scale': 1/DPS/ABM_Vol, 'ylim': (), 'type': 'inc', 'col':'prolif.cancerCell', 'title': 'Cancer cell prolif, '+r'$1/(day\cdot mm^{3})$'},
                {'loc': (0, 4), 'scale': 1/DPS/ABM_Vol, 'ylim': (), 'type': 'inc', 'col':'death.cancerCell', 'title': 'Cancer cell death, '+r'$1/(day\cdot mm^{3})$'},
                {'loc': (0, 5), 'scale': 1/DPS/ABM_Vol, 'ylim': (),  'type': 'inc', 'col':'killed_by_t', 'title': 'Cancer cell killed, '+r'$1/(day\cdot mm^{3})$'},
                {'loc': (1, 0), 'scale': 1/ABM_Vol, 'ylim': (), 'type': 'raw', 'col':'agentCount.CD8.effector', 'title': 'CD8+ T cell, eff, '+r'$1/mm^{3}$'},
                {'loc': (1, 1), 'scale': 1/ABM_Vol, 'ylim': (), 'type': 'raw', 'col':'agentCount.CD8.cytotoxic', 'title': 'CD8+ T cell, cyt, '+r'$1/mm^{3}$'},
                {'loc': (1, 2), 'scale': 1/ABM_Vol, 'ylim': (), 'type': 'raw', 'col':'agentCount.CD8.suppressed', 'title': 'CD8+ T cell, exh, '+r'$1/mm^{3}$'},
                {'loc': (1, 3), 'scale': 1/DPS/ABM_Vol, 'ylim': (), 'type': 'inc', 'col':'recruit.CD8.effector', 'title': 'CD8+ T cell, recruitment, '+r'$1/(day\cdot mm^{3})$'},
                {'loc': (1, 4), 'scale': 1/ABM_Vol, 'ylim': (), 'type': 'raw', 'col':'agentCount.Treg.default', 'title': 'Treg, '+r'$1/mm^{3}$'},
                {'loc': (1, 5), 'scale': 1, 'ylim': (0, 1.1), 'type': 'raw', 'col':'H_PD1_PDL1', 'title': 'H_PD1_PDL1'},
               ]
    
    # figure
    #fig, axs = plt.subplots(5,5, constrained_layout=True)
    num_row = 2*num_window
    num_col = 6
    fig, axs = plt.subplots(num_row, num_col)
    
    for s in samples:
        sim_id = get_sim_id(db, group, s, treat, 1) # None if no such simulation archived
        if sim_id is not None:
            
            #K_C1_PDLX_Teff_PD1 = get_param_val(db, 'QSP/init_value/Parameter/K_C1_PDLX_Teff_PD1',  params, group, s, treat) # item
            #n_PD1_PDLX = get_param_val(db, 'QSP/init_value/Parameter/n_PD1_PDLX',  params, group, s, treat) # dimensionless
            for j, window in enumerate(windows):
                for i, var in enumerate(var_list):
                    ax_row, ax_col = var['loc']
                    ax = axs[j*2+ax_row, ax_col]
                    #print(var['col'])
                    if var['type'] == 'raw':
                        y = get_abm_time_series(db, group, s, treat, 1, window, var['col'])
                    elif var['type'] == 'inc':
                        y = get_abm_time_series_inc(db, group, s, treat, 1, window, var['col'])
                    ax.plot(t_qsp_day, y*var['scale'],  alpha=alpha, color = colors[s-1])
                    ax.set_title(var['title'])
                    ax.set_xlabel('t, day')
                    ax.set_yscale('log')
                    if var['ylim'] and not use_auto_ylim:
                        ax.set_ylim(*var['ylim'])
                        ax.set_yscale('linear')
                    if var['loc'][1] == 0:
                        ax.set_ylabel(window)
            
    for r in range(num_col):
        axs[num_row-1,0].set_xlabel('t, day')
    
    fig.set_size_inches(24,16)
    fig.tight_layout()
    return fig

#%% summary species with replications
def get_all_rep(db, gid, tid, sid, reps, vartype, var, window = None):
    VAR_QSP = 'qsp'
    VAR_ABM_RAW = 'abm'
    VAR_ABM_INC = 'abm_inc'
    
    sim_id_0 = db.fetch(qdb.TABLE_NAME_QSP, ['UID'], '')[0][0]
    n_time = len(get_qsp_time_series_by_sim_id(db, sim_id_0, ['time']).flatten())
    res_reps = np.zeros((len(reps), n_time))
    n = 0
    
    for rid in reps:
        sim_id = get_sim_id(db, gid, sid, tid, rid)
        if sim_id is not None: # result exist
            res = None
            if vartype == VAR_QSP:
                res = get_qsp_time_series(db, gid, sid, tid, rid, var)
            elif vartype == VAR_ABM_RAW:
                res = get_abm_time_series(db, gid, sid, tid, rid, window, var)
            elif vartype == VAR_ABM_INC:
                res = get_abm_time_series_inc(db, gid, sid, tid, rid, window, var)
            else:
                raise ValueError('varType not recognized: {}'.format(vartype))
            res_reps[n, :] = res
            n += 1
    return res_reps[:n, :]

def get_diam_rep(db, gid, sid, tid, reps, t_nivo_start, weight_qsp):
    res_reps = None
    vol_reps = None
    n = 0
    for rid in reps:
        sim_id = get_sim_id(db, gid, sid, tid, rid)
        if sim_id is not None: # result exist
            d_chg_time, vol = get_tum_size(db, gid, sid, tid, rid, t_nivo_start, weight_qsp)
            if res_reps is None:
                res_reps = np.zeros((len(reps), len(d_chg_time)))
                vol_reps = np.zeros((len(reps), len(d_chg_time)))
            res_reps[n, :] = d_chg_time
            vol_reps[n, :] = vol
            n += 1            
    return res_reps[:n, :], np.mean(vol_reps, axis=0)

def plot_rep(ax, t, res, color, alpha):
    sem = np.std(res, axis=0, ddof=1) / np.sqrt(res.shape[0])
    y = np.mean(res, axis=0)
    ax.fill_between(t, y+sem, y-sem, alpha=alpha, color = color, linewidth=0.0)
    ax.plot(t, y, color = color)
    return
 
def plot_avg_summary(db, params, groups, treat, sid, reps, alpha):
    
    #colors = ['tab:red', 'tab:green', 'tab:blue', 'tab:olive']
    colors = ['tab:red', 'olivedrab', 'tab:blue', 'gold']
    
    NIVO_START_DAY = get_param_val(db, 'QSP/init_value/Parameter/t_init_Nivo', params, 1, 1, 1)
    SLICE_PER_DAY = SEC_PER_DAY / get_param_val(db, 'ABM/Environment/SecPerSlice', params, 1, 1, 1)
    t_nivo_start = NIVO_START_DAY*SLICE_PER_DAY
    
    sim_id_0 = db.fetch(qdb.TABLE_NAME_QSP, ['UID'], '')[0][0]
    t = get_qsp_time_series_by_sim_id(db, sim_id_0, ['time']).flatten()/SLICE_PER_DAY
    
    V2A = 0.005+0.01-2*0.0025
    V_C = get_param_val(db, 'QSP/init_value/Parameter/vol_cent', params, 1, 1, 1) * 1e6 #mm^3
    DPS = t[1]-t[0]
    
    XSIZE = get_param_val(db, 'ABM/Environment/Tumor/XSize', params, 1, 1, 1)
    YSIZE = get_param_val(db, 'ABM/Environment/Tumor/YSize', params, 1, 1, 1)
    ZSIZE = get_param_val(db, 'ABM/Environment/Tumor/ZSize', params, 1, 1, 1)
    VoxelSize = get_param_val(db, 'ABM/Environment/Tumor/VoxelSize', params, 1, 1, 1)
    V_ABM = XSIZE * YSIZE * ZSIZE * VoxelSize**3 * 1e-9
    
    AVO = get_param_val(db, 'QSP/init_value/Parameter/AvogadroN', params, 1, 1, 1)
    
    fig, axs = plt.subplots(4,6)
    config = [# row 1: QSP
                {'title':'Tumor diameter\n% change', 'scale':'linear', 'ylim': (-.5, 2)},
                {'title':'Teff, Tumor\n$1/mm^{2}$', 'scale':'linear', 'ylim': ()},
                {'title':'Treg, Tumor\n$1/mm^{2}$', 'scale':'linear', 'ylim': ()},
                {'title':'Teff, Blood\n$1/mm^{3}$', 'scale':'linear', 'ylim': ()},
                {'title':'Treg, Blood\n$1/mm^{3}$', 'scale':'linear', 'ylim': ()},
                {'title':'Nivo, Tumor\nnM', 'scale':'linear', 'ylim': ()},
                # row 2: ABM
                {'title':'Tcyt, core\n$1/mm^{3}$', 'scale':'log', 'ylim': (1e1, 5e4)},
                {'title':'Treg, core\n$1/mm^{3}$', 'scale':'log', 'ylim': (1e1, 5e4)},
                {'title':'CC killed, core\n$1/(mm^{3}\cdot day)$', 'scale':'log', 'ylim': (10, 2e4)},
                {'title':'Tcyt, IF\n$1/mm^{3}$', 'scale':'log', 'ylim': (1e1, 5e4)},
                {'title':'Treg, IF\n$1/mm^{3}$', 'scale':'log', 'ylim': (1e1, 5e4)},
                {'title':'CC killed, IF\n$1/(mm^{3}\cdot day)$', 'scale':'log', 'ylim': (10,1e4)},
                # row 3: cytotoxicity
                {'title':'PDL1 inhibition', 'scale':'linear', 'ylim': ()},
                {'title':'CC fractional change\n 1/day', 'scale':'linear', 'ylim': ()},
                {'title':'Fraction CC killed\n 1/day', 'scale':'linear', 'ylim': ()},
                {'title':'Teff/Total ratio, %', 'scale':'linear', 'ylim': ()},
                {'title':'Treg/Teff ratio, %', 'scale':'linear', 'ylim': ()},
                {'title':'Teff frational exhaustion\nby Treg, 1/day', 'scale':'linear', 'ylim': ()},
                # row 4: immune response
                {'title':'Tumor antigen\npM', 'scale':'linear', 'ylim': ()},
                {'title':'LN mAPC', 'scale':'linear', 'ylim': ()},
                {'title':'pMHC, mAPC', 'scale':'linear', 'ylim': ()},
                {'title':'CD8 activated', 'scale':'linear', 'ylim': ()},
                {'title':'LN IL2, nM', 'scale':'linear', 'ylim': ()},
                {'title':'CD8 division\nIL2 induced', 'scale':'linear', 'ylim': ()}
        ]
    
    for i, gid in enumerate(groups):
        c = colors[i]
        ## QSP
        weight_qsp = get_param_val(db, 'QSP/simulation/weight_qsp',  params, gid, 1, treat)
        Ccell = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.C1')
        Ccell[Ccell<1]=1
        T1_tum = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.Teff_1_0')
        Tr_tum = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.Treg')
        Te_tum = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.Teff_exhausted')
        
        # tumor diameter(spider)
        ax = axs[0,0]
        res, vol_time = get_diam_rep(db, gid, sid, treat, reps, t_nivo_start, weight_qsp)
        plot_rep(ax, t-NIVO_START_DAY, res, color=c, alpha=alpha)
        ax.set_xlim(0, t[-1]-NIVO_START_DAY)
            
        # Tumor Teff 
        ax = axs[0,1]
        plot_rep(ax, t, (T1_tum+Te_tum)/vol_time*V2A, color=c, alpha=alpha)
      
        # tumor Treg
        ax = axs[0,2]
        plot_rep(ax, t, (Tr_tum)/vol_time*V2A, color=c, alpha=alpha)
        
        # blood Teff
        ax = axs[0,3]
        Te_Cent = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Cent.Teff_1_0')/V_C
        plot_rep(ax, t, Te_Cent, color=c, alpha=alpha)
        
        # blood Treg
        ax = axs[0,4]
        Tr_Cent = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Cent.Treg')/V_C
        plot_rep(ax, t, Tr_Cent, color=c, alpha=alpha)
        
        # Nivo
        ax = axs[0,5]
        res = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.Nivo') * 1e9
        plot_rep(ax, t, res, color=c, alpha=alpha)
        
        ## ABM
        # Teff, Core
        ax = axs[1,0]
        res = get_all_rep(db, gid, treat, sid, reps, 'abm', 'agentCount.CD8', window = 'core')/V_ABM
        plot_rep(ax, t, res, color=c, alpha=alpha)
        # Treg, Core
        ax = axs[1,1]
        res = get_all_rep(db, gid, treat, sid, reps, 'abm', 'agentCount.Treg', window = 'core')/V_ABM
        plot_rep(ax, t, res, color=c, alpha=alpha)
        # Cancer killed, Core
        ax = axs[1,2]
        res = get_all_rep(db, gid, treat, sid, reps, 'abm_inc', 'killed_by_t', window = 'core')/DPS/V_ABM
        plot_rep(ax, t, res, color=c, alpha=alpha)

        # Teff, Core
        ax = axs[1,3]
        res = get_all_rep(db, gid, treat, sid, reps, 'abm', 'agentCount.CD8', window = 'margin')/V_ABM
        plot_rep(ax, t, res, color=c, alpha=alpha)
        # Treg, Core
        ax = axs[1,4]
        res = get_all_rep(db, gid, treat, sid, reps, 'abm', 'agentCount.Treg', window = 'margin')/V_ABM
        plot_rep(ax, t, res, color=c, alpha=alpha)
        # Cancer killed, Core
        ax = axs[1,5]
        res = get_all_rep(db, gid, treat, sid, reps, 'abm_inc', 'killed_by_t', window = 'margin')/DPS/V_ABM
        plot_rep(ax, t, res, color=c, alpha=alpha)

        # Cytotoxicity        
        k_C_death_by_T = get_param_val(db, 'QSP/init_value/Parameter/k_C_death_by_T',  params, gid, sid, treat) # 1/day
        k_Teff_death_by_C = get_param_val(db, 'QSP/init_value/Parameter/k_Teff_death_by_C',  params, gid, sid, treat)# 1/day
        k_Teff_inhibBy_Treg = get_param_val(db, 'QSP/init_value/Parameter/k_Teff_inhibBy_Treg',  params, gid, sid, treat) #1/day
        K_C1_PDLX_Teff_PD1 = get_param_val(db, 'QSP/init_value/Parameter/K_C1_PDLX_Teff_PD1',  params, gid, sid, treat) # item
        n_PD1_PDLX = get_param_val(db, 'QSP/init_value/Parameter/n_PD1_PDLX',  params, gid, sid, treat) # dimensionless
        k_C_death = get_param_val(db, 'QSP/init_value/Parameter/k_C_death',  params, gid, sid, treat) # 1/day
        k_C_growth = get_param_val(db, 'QSP/init_value/Parameter/k_C_growth',  params, gid, sid, treat) # 1/day
        vol_tum_max = get_param_val(db, 'QSP/init_value/Parameter/vol_tum_max',  params, gid, sid, treat) # L
        D_C = get_param_val(db, 'QSP/init_value/Parameter/D_C',  params, gid, sid, treat) # micron
        f_vol_tum = get_param_val(db, 'QSP/init_value/Parameter/f_vol_tum',  params, gid, sid, treat) # L
        
        vol_C = 3.1416/6*D_C**3 * 1e-15 # L
        K_C_MAX = vol_tum_max/vol_C*(1-f_vol_tum)
        
        tum_C = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.C1')
        tum_Treg = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.Treg')
        tum_Teff = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.Teff_1_0')
        tum_T_exh = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.Teff_exhausted')
        tum_all = tum_C + tum_Teff + tum_T_exh + tum_Treg # exhausted T cell accounted
        PDLX_PD1 = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.C1_PDL1_Teff_PD1') + \
            get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.C1_PDL2_Teff_PD1')
        
        
        H_PDLX_PD1_n = f_hill(PDLX_PD1, K_C1_PDLX_Teff_PD1, n_PD1_PDLX)
        H_PDLX_PD1_1 = f_hill(PDLX_PD1, K_C1_PDLX_Teff_PD1, 1)
        
        dcdt_g = k_C_growth * (1-tum_C/K_C_MAX)
        dcdt_d = k_C_death
        dcdt_k = k_C_death_by_T * tum_Teff/tum_all *(1-H_PDLX_PD1_n)
            
        # PDL1 inhibition
        ax = axs[2,0]
        plot_rep(ax, t, H_PDLX_PD1_n, color=c, alpha=alpha)
        # Fractional appaarent growth
        ax = axs[2,1]
        plot_rep(ax, t, dcdt_g - dcdt_d - dcdt_k, color=c, alpha=alpha)
        # Fractional death by killing
        ax = axs[2,2]
        plot_rep(ax, t, dcdt_k, color=c, alpha=alpha)
        # Teff/Total
        ax = axs[2,3]
        plot_rep(ax, t, 100 * tum_Teff/tum_all, color=c, alpha=alpha)
        # Treg/Teff ratio        
        ax = axs[2,4]
        plot_rep(ax, t, 100 * Tr_tum/(T1_tum+Te_tum), color=c, alpha=alpha)
        # Exhaustion by Treg
        ax = axs[2,5]
        inhib_by_reg = k_Teff_inhibBy_Treg *tum_Treg/tum_all*(1+H_PDLX_PD1_n)
        plot_rep(ax, t, inhib_by_reg, color=c, alpha=alpha)
        
        # Immune
        n_aT_prolif_IL20 = get_param_val(db, 'QSP/init_value/Parameter/n_aT_prolif_IL20',  params, gid, sid, treat) # dimensionless
        f_vol_LN =get_param_val(db, 'QSP/init_value/Parameter/f_vol_LN',  params, gid, sid, treat) # dimensionless
        K_IL2_Teff = get_param_val(db, 'QSP/init_value/Parameter/K_IL2_Teff',  params, gid, sid, treat) # nmole/L
        
        IL2 = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'LN.IL2')
        n_aT_prolif_IL2 = n_aT_prolif_IL20 * (IL2 / f_vol_LN) / (IL2 / f_vol_LN + K_IL2_Teff) + 0.001
        
        # Tumor Ag
        ax = axs[3,0]
        D1_0 = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'Tum.D1_0')*1e12
        plot_rep(ax, t, D1_0, color=c, alpha=alpha)
        # LN mAPC
        ax = axs[3,1]
        mAPC = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'LN.mAPC')
        plot_rep(ax, t, mAPC, color=c, alpha=alpha)
        # mAPC pMHC
        ax = axs[3,2]
        pMHC = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'mAPC_Surf.p1_0_M1')*AVO
        plot_rep(ax, t,  pMHC, color=c, alpha=alpha)
        # CD8 activated
        ax = axs[3,3]
        aT = get_all_rep(db, gid, treat, sid, reps, 'qsp', 'LN.aT_1_0')
        plot_rep(ax, t,  aT, color=c, alpha=alpha)
        # LN IL2    
        ax = axs[3,4]
        plot_rep(ax, t, IL2, color=c, alpha=alpha)
        # number of division
        ax = axs[3,5]
        plot_rep(ax, t,  n_aT_prolif_IL2, color=c, alpha=alpha)
    
    from matplotlib.ticker import FormatStrFormatter    
    for i, ax in enumerate(axs.flatten()):
        if config[i]['ylim']:
            ax.set_ylim(*config[i]['ylim'])
        ax.set_yscale(config[i]['scale'])
        if 'log' != config[i]['scale']:
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
        ax.set_title(config[i]['title'])
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    return fig    

#%% Other supporting functions
# hill equation
def f_hill(L, k, n):
    return L**n/(L**n + k**n)