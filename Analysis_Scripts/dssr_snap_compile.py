#! /usr/bin/python
import os
import numpy as np
import pandas as pd

class Snap:
    def __init__(self, name):
        self.name = name
    
    def extract_data(self, filename):
        '''
        Read dssr-snap output into list
        '''
        infile=open(filename, "r")
        data_list=infile.readlines()
        infile.close()
        data_list=[i.rstrip('\n') for i in data_list]
        return data_list
    
    def nucleotide_aa_dataframe(self, data_list):
        '''
        Go through dssr-snap list and extract the nucleotide/amino-acid interactions section.
        Put section into formatted dataframe
        '''
        for line in data_list:
            if 'nucleotide/amino-acid interactions' in line:
                A = data_list.index(line)
                B = int(line.split()[2])
                lst=[i.split() for i in data_list[A+1:A+1+B+1]]
                ntaa = pd.DataFrame(lst[1:], columns=['index']+lst[0]).reset_index(drop=True).set_index('index')
                ntaa[['Tdst','Rdst','Tx','Ty','Tz','Rx','Ry','Rz']] = ntaa[['Tdst','Rdst','Tx','Ty','Tz','Rx','Ry','Rz']].astype(float)
                del A, B, lst
        return ntaa
    
    def basepair_aa_dataframe(self, data_list):
        '''
        Go through dssr-snap list and extract the base-pair/amino-acid interactions section.
        Put section into formatted dataframe
        '''
        for line in data_list:
            if 'base-pair/amino-acid interactions' in line:
                A = data_list.index(line)
                B = int(line.split()[2])
                lst=[i.split() for i in data_list[A+1:A+1+B+1]]
                bpaa = pd.DataFrame(lst[1:], columns=['index']+lst[0]).reset_index(drop=True).set_index('index')
                bpaa[['Tdst','Rdst','Tx','Ty','Tz','Rx','Ry','Rz']] = bpaa[['Tdst','Rdst','Tx','Ty','Tz','Rx','Ry','Rz']].astype(float)
                del A, B, lst
        return bpaa
    
    def phosphate_aa_hbond_dataframe(self, data_list):
        '''
        Go through dssr-snap list and extract the phosphate/amino-acid interactions section.
        Put section into formatted dataframe
        '''
        for line in data_list:
            if 'phosphate/amino-acid H-bonds' in line:
                A = data_list.index(line)
                B = int(line.split()[2])
                lst=[i.split() for i in data_list[A+1:A+1+B+1]]
                hpaa = pd.DataFrame(lst[1:], columns=['index']+lst[0]).reset_index(drop=True).set_index('index')
                hpaa['dist']=hpaa['dist'].astype(float)
                del A, B, lst
        return hpaa
    
    def basepair_aa_hbond_dataframe(self, data_list):
        '''
        Go through dssr-snap list and extract the base/amino-acid H-bonds section.
        Put section into formatted dataframe
        '''
        for line in data_list:
            if 'base/amino-acid H-bonds' in line:
                A = data_list.index(line)
                B = int(line.split()[2])
                lst=[i.split() for i in data_list[A+1:A+1+B+1]]
                hbpp = pd.DataFrame(lst[1:], columns=['index']+lst[0]).reset_index(drop=True).set_index('index')
                hbpp['dist']=hbpp['dist'].astype(float)
                del A, B, lst
        return hbpp
    
    def base_aa_dataframe(self, data_list):
        '''
        Go through dssr-snap list and extract the base/amino-acid pairs section.
        Put section into formatted dataframe
        '''
        for line in data_list:
            if 'base/amino-acid pairs' in line:
                A = data_list.index(line)
                B = int(line.split()[2])
                lst=[i.split() for i in data_list[A+1:A+1+B+1]]
                baa = pd.DataFrame(lst[1:], columns=['index']+lst[0]).reset_index(drop=True).set_index('index')
                baa[['vertical-distance','plane-angle']]=baa[['vertical-distance','plane-angle']].astype(float)
                del A, B, lst
        return baa
    
    def sugar_aa_hbond_dataframe(self, data_list):
        '''
        Go through dssr-snap list and extract the base/amino-acid pairs section.
        Put section into formatted dataframe
        '''
        for line in data_list:
            if 'sugar/amino-acid H-bonds' in line:
                A = data_list.index(line)
                B = int(line.split()[2])
                lst=[i.split() for i in data_list[A+1:A+1+B+1]]
                hsaa = pd.DataFrame(lst[1:], columns=['index']+lst[0]).reset_index(drop=True).set_index('index')
                hsaa['dist']=hsaa['dist'].astype(float)
                del A, B, lst
        return hsaa

