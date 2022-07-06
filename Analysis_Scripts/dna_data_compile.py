#! /usr/bin/python
import os, sys
import math, scipy
import numpy as np
import pandas as pd
path = os.getcwd()
if "Young_Research" in path:
    initialpath = "C:\\Users\\Young_Research\\Documents\\Rutgers\\Research\\InitialConditions"
    ffpath = "C:\\Users\\Young_Research\\Documents\\Rutgers\\Research\\ForceFields"
    sys.path.append("C:\\Users\\Young_Research\\Documents\\Rutgers\\Research\\Scripts\\analysis\\")
else:
    initialpath= "/home/rty10/Documents/InitialConditions"
    ffpath = "/home/rty10/Documents/ForceFields"
    sys.path.append("/home/rty10/Documents/Scripts/analysis/")
import dna_analysis

ematrix_dict = {0:'tilt', 1:'roll', 2:'twist', 3:'shift', 4:'slide', 5:'rise'}


def file_to_list(filename):
    infile = open(filename, 'r')
    data = infile.readlines()
    data = [i.rstrip('\n') for i in data]
    infile.close()
    return data
    

def log_file_data(filename, df_index, dataframe):
    '''
    Collect initial energy, final energy, and the energy matrix data
    '''
    ematrix_dict = {0:'tilt', 1:'roll', 2:'twist', 3:'shift', 4:'slide', 5:'rise'}
    log_file_list = file_to_list(filename)
    
    for i in range(0, len(log_file_list)):
        if 'initial energy:' in log_file_list[i]:
            dataframe.at[df_index, 'eo'] = float( log_file_list[i].split(':')[1] )
        elif 'final energy:' in log_file_list[i]:
            dataframe.at[df_index, 'eopt'] = float( log_file_list[i].split(':')[1] )
    
    ematrix = log_file_list[-7:-1]
    
    for i in range(0, len(ematrix)):
        ematrix[i] = ematrix[i].replace('{', '').replace('}', '').split(',')
        for j in range(0, len(ematrix[i])):
            ematrix[i][j] = float(ematrix[i][j])
        dataframe.at[df_index, 'eopt_'+ematrix_dict[i]] = ematrix[i][i]
    
    del log_file_list, ematrix_dict
    return dataframe


def topology_file_data(filename, df_index, dataframe):
    '''
    '''
    topo_file_list = file_to_list(filename)
    topo_file_list = topo_file_list[-4:]
    for i in range(0, len(topo_file_list)):
        if 'Wr ' in topo_file_list[i]:
            dataframe.at[df_index, 'wr'] = float( topo_file_list[i].split('=')[1] )
        elif 'Tw ' in topo_file_list[i]:
            dataframe.at[df_index, 'tw'] = float( topo_file_list[i].split('=')[1] )
        elif 'Lk ' in topo_file_list[i]:
            dataframe.at[df_index, 'lk'] = int( topo_file_list[i].split('=')[1] )
    del topo_file_list[:]
    return dataframe

def refframe_file_data(filename, df_index, dataframe, LENGTH_TOTAL, LENGTH_LOOP, start_bp_id, midloop_bp_id, end_bp_id):
    '''
    '''
    rfdf    = dna_analysis.df_read_refframe_origins(filename)
    rfnorms = dna_analysis.df_read_refframe_normals(filename)
    
    if len(rfdf) != LENGTH_TOTAL:
        rfdf = rfdf[:-(abs(len(rfdf)-LENGTH_TOTAL))]
    if len(rfnorms) != LENGTH_TOTAL:
        rfnorms = rfnorms[:-(abs(len(rfnorms)-LENGTH_TOTAL))]
    
    dataframe.at[df_index, 'rg'] = dna_analysis._gyration_origindf(rfdf)
    
    dataframe.at[df_index, 'enex_dis']   = dna_analysis._vec_distance(rfdf.loc[end_bp_id] - rfdf.loc[start_bp_id])
    dataframe.at[df_index, 'enex_angle'] = dna_analysis._loop_origins_omega_angles(rfnorms.loc[end_bp_id], rfnorms.loc[start_bp_id])
	
    WE, WA, WB, WG = dna_analysis._loop_westcott_angles(rfdf, rfnorms, LENGTH_LOOP, start_bp_id, midloop_bp_id, end_bp_id)
    dataframe.at[df_index, 'loop_eta']=WE
    dataframe.at[df_index, 'loop_alpha'] = WA
    dataframe.at[df_index, 'loop_beta']=WB
    dataframe.at[df_index, 'loop_gamma']=WG
    del WE, WA, WB, WG
    
    P1, P2, P3 = dna_analysis._loop_pinching(rfdf, LENGTH_LOOP, start_bp_id, end_bp_id)
    dataframe.at[df_index, 'loop_pinchpts']=P1
    dataframe.at[df_index, 'loop_pinchdis']=P2
    dataframe.at[df_index, 'loop_pinchlen']=P3
    del P1, P2, P3
    
    del rfdf, rfnorms
    return dataframe
