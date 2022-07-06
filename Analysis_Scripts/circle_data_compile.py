#! /usr/bin/python
import os, sys
import math, scipy
import numpy as np
import pandas as pd
initialpath = "C:\\Users\\Young_Research\\Documents\\Rutgers\\Research\\InitialConditions"
ffpath = "C:\\Users\\Young_Research\\Documents\\Rutgers\\Research\\ForceFields"
sys.path.append("C:\\Users\\Young_Research\\Documents\\Rutgers\\Research\\Scripts\\analysis\\")
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
        elif 'elastic ->' in log_file_list[i]:
            dataframe.at[df_index, 'eopt_elastic'] = float( log_file_list[i].split('->')[1] )
        elif 'electrostatic ->' in log_file_list[i]:
            dataframe.at[df_index, 'eopt_electrostatic'] = float( log_file_list[i].split('->')[1] )
    
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

def refframe_file_data(filename, df_index, dataframe, Ncirc):
    '''
    '''
    rfdf    = dna_analysis.df_read_refframe_origins(filename)
    rfnorms = dna_analysis.df_read_refframe_normals(filename)
    
    
    if len(rfdf) != Ncirc:
        rfdf = rfdf[:-(abs(len(rfdf)-Ncirc))]
    if len(rfnorms) != Ncirc:
        rfnorms = rfnorms[:-(abs(len(rfnorms)-Ncirc))]
    
    dataframe.at[df_index, 'rg'] = dna_analysis._gyration_origindf(rfdf)
    
    pca = dna_analysis._pca_origindf(rfdf)
    dataframe.at[df_index, 'ellipticity'] = dna_analysis._ellipticity(pca)
    del pca
    
    del rfdf, rfnorms
    return dataframe
