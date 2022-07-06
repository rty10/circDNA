#! /usr/bin/python
import os, sys, subprocess, shutil
import scipy
from scipy import spatial
import numpy as np
import pandas as pd

path = os.getcwd()

def dna_data(refframe_file, pdb_file):
    """
    Function that will read both a reference frame file and a file of phosphate atoms from a pdb file.
    Goal will be to return a pandas DataFrame with information such as:
    - base-pair origin
    - base-pair's direction of its minor groove
    - the phosphate positions on both its coding (Pi) and complementary (pi) strands
    """
    infile1   = open(refframe_file, 'r')
    infile2   = open(pdb_file, 'r')
    rfdata    = infile1.readlines()
    PDBdata   = infile2.readlines()
    infile1.close()
    infile2.close()
    
    df = pd.DataFrame(columns=["Ox","Oy","Oz",
                               "Tox","Toy","Toz",
                               "-d1x","-d1y","-d1z",
                               "px","py","pz",
                               "qx","qy","qz"])    
    Ncirc = int(rfdata[0].split()[0])
    rfdata  = [i.split() for i in rfdata]
    Odata = []
    Xdata = []
    for i in range(0, Ncirc):
        Odata.append(rfdata[(5*i)+2])
        Xdata.append([rfdata[(5*i)+3][0],
                      rfdata[(5*i)+4][0],
                      rfdata[(5*i)+5][0]])
    for i in range(0, len(Odata)):
        df.loc[i, ["Ox","Oy","Oz"]]     = float(Odata[i][0]), float(Odata[i][1]), float(Odata[i][2])
        df.loc[i, ["-d1x","-d1y","-d1z"]] = -1*float(Xdata[i][0]), -1*float(Xdata[i][1]), -1*float(Xdata[i][2])
    
    com = np.array([float((df['Ox'].sum())/Ncirc),
                    float((df['Oy'].sum())/Ncirc),
                    float((df['Oz'].sum())/Ncirc)])
    for i in range(0, len(Odata)):
        df.loc[i, ["Tox","Toy","Toz"]] = (np.array(df.loc[i, ["Ox","Oy","Oz"]]) - com)
    
    Pdata = []
    for i in range(0, len(PDBdata)):
        if " P " in PDBdata[i]:
            Pdata.append(PDBdata[i])
    Pdata = [[i[30:38],i[38:46],i[46:54]] for i in Pdata]
    # if this opt as a circle, 1 and Ncirc+1 will be identical    
    if len(Pdata) == 2*Ncirc + 2:
        for i in range(0, len(Odata)):
            df.loc[i, ["px","py","pz"]] = float(Pdata[i][0]), float(Pdata[i][1]), float(Pdata[i][2])
            j = ((2*Ncirc + 2)-i)-1
            df.loc[i, ["qx","qy","qz"]] = float(Pdata[j][0]), float(Pdata[j][1]), float(Pdata[j][2])
    del rfdata, PDBdata, Odata, Xdata, Pdata
    df = circ_minor_groove(df)
    df = circ_major_groove(df)
    return df

def circ_minor_groove(dataframe):
    """
    To calculate the width of the minor groove for a circular structure.
    For step 'k', get the average of distances between a pair of phosphates offset by m=-3.
    A: distance from coding p(k+1) to complementary q(k-2)
    B: distance from coding p(k+2) to complementary q(k-1)
    *** For this function, the p phosphates will be complementary to the k-1th step's bp number
    *** for k=10 and circular N=150, p((k+1)+2)=p(k+3)=p([13])
    """
    Nseq = len(dataframe)
    for k in range(0, Nseq):
        if k == 0:
            Pk1 = np.array([dataframe.loc[2,      ["px","py","pz"]]])
            Pk2 = np.array([dataframe.loc[3,      ["px","py","pz"]]])
            pk1 = np.array([dataframe.loc[Nseq-1, ["qx","qy","qz"]]])
            pk2 = np.array([dataframe.loc[Nseq-2, ["qx","qy","qz"]]])
        elif k == 1:
            Pk1 = np.array([dataframe.loc[3,      ["px","py","pz"]]])
            Pk2 = np.array([dataframe.loc[4,      ["px","py","pz"]]])
            pk1 = np.array([dataframe.loc[0,      ["qx","qy","qz"]]])
            pk2 = np.array([dataframe.loc[Nseq-1, ["qx","qy","qz"]]])
        elif k == Nseq-3:
            Pk1 = np.array([dataframe.loc[k+2,    ["px","py","pz"]]])
            Pk2 = np.array([dataframe.loc[0,      ["px","py","pz"]]])
            pk1 = np.array([dataframe.loc[k-1,    ["qx","qy","qz"]]])
            pk2 = np.array([dataframe.loc[k-2,    ["qx","qy","qz"]]])
        elif k == Nseq-2:
            Pk1 = np.array([dataframe.loc[0,      ["px","py","pz"]]])
            Pk2 = np.array([dataframe.loc[1,      ["px","py","pz"]]])
            pk1 = np.array([dataframe.loc[k-1,    ["qx","qy","qz"]]])
            pk2 = np.array([dataframe.loc[k-2,    ["qx","qy","qz"]]])
        elif k == Nseq-1:
            Pk1 = np.array([dataframe.loc[1,      ["px","py","pz"]]])
            Pk2 = np.array([dataframe.loc[2,      ["px","py","pz"]]])
            pk1 = np.array([dataframe.loc[k-1,    ["qx","qy","qz"]]])
            pk2 = np.array([dataframe.loc[k-2,    ["qx","qy","qz"]]])
        else:
            Pk1 = np.array([dataframe.loc[k+2,    ["px","py","pz"]]])
            Pk2 = np.array([dataframe.loc[k+3,    ["px","py","pz"]]])
            pk1 = np.array([dataframe.loc[k-1,    ["qx","qy","qz"]]])
            pk2 = np.array([dataframe.loc[k-2,    ["qx","qy","qz"]]])
        
        A = sp.spatial.distance.euclidean(Pk1, pk2)
        B = sp.spatial.distance.euclidean(Pk2, pk1)
        mg = (1/2)*( A + B )
        dataframe.loc[k, "W-min"] = mg
    del Nseq
    return dataframe

def circ_major_groove(dataframe):
    """
    To calculate the width of the major groove for a circular structure.
    For step 'k', get the average of distances between a pair of phosphates offset by m=4.
    mg = distance from coding p((k+1)-2)=p(k-1) to complementary q(k+2)
    """
    Nseq = len(dataframe)
    for k in range(0, Nseq):
        if k == 0:
            Pk2 = np.array([dataframe.loc[Nseq-1, ["px","py","pz"]]])
            pk2 = np.array([dataframe.loc[2,      ["qx","qy","qz"]]])
        elif k == 1:
            Pk2 = np.array([dataframe.loc[0,      ["px","py","pz"]]])
            pk2 = np.array([dataframe.loc[3,      ["qx","qy","qz"]]])
        elif k == Nseq-2:
            Pk2 = np.array([dataframe.loc[k-1,    ["px","py","pz"]]])
            pk2 = np.array([dataframe.loc[0,      ["qx","qy","qz"]]])
        elif k == Nseq-1:
            Pk2 = np.array([dataframe.loc[k-1,    ["px","py","pz"]]])
            pk2 = np.array([dataframe.loc[1,      ["qx","qy","qz"]]])
        else:
            Pk2 = np.array([dataframe.loc[k-1,    ["px","py","pz"]]])
            pk2 = np.array([dataframe.loc[k+2,    ["qx","qy","qz"]]])
            
        mg = sp.spatial.distance.euclidean(Pk2, pk2)
        dataframe.loc[k, "W-maj"] = mg
    del Nseq
    return dataframe

