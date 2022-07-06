"""
    Stefjord Todolli
    July 12, 2020
"""

import os, sys, ast
import pathlib, shutil, tempfile
import numpy as np
import pandas as pd
import subprocess as sub
from typing import List

def df_read_refframe_origins(refframe_filename):
    '''
    '''
    df   = pd.DataFrame(columns=['x','y','z'])
    infile = open(refframe_filename, 'r')
    rfdata = infile.readlines()
    infile.close()
    rfdata = [i.rstrip('\n').replace('  # origin','').replace('  # x-axis','').replace('  # y-axis','').replace('  # z-axis','').split() for i in rfdata]
    N = int(rfdata[0][0])
    for j in range(0, N):
        df.at[j, ['x','y','z']]  = rfdata[ 5*j + 2]
    df = df.astype('float')
    df.index += 1
    return df

def df_read_refframe_vector(refframe_filename, vector):
    '''
    '''
    vector_dict={'x':3, 'y':4, 'z':5}
    df   = pd.DataFrame(columns=['x','y','z'])
    infile = open(refframe_filename, 'r')
    rfdata = infile.readlines()
    infile.close()
    rfdata = [i.rstrip('\n').replace('  # origin','').replace('  # x-axis','').replace('  # y-axis','').replace('  # z-axis','').split() for i in rfdata]
    N = int(rfdata[0][0])
    for j in range(0, N):
        df.at[j, ['x','y','z']]  = rfdata[ 5*j + vector_dict[vector]]
    df = df.astype('float')
    df.index += 1
    return df


class GeoToC:

    def __init__(self, executable: str = "/usr/local/bin/GeoToC_CLU"):
        self.exe = executable

    def analyze_twdensity_curves(self, curves: List) -> List[float]:
        originaldir = pathlib.Path.cwd()
        tempdir = tempfile.mkdtemp()
        os.chdir(tempdir)
        np.savetxt("curve1.tsv", np.array([c[0] for c in curves]), fmt=['%.5f', '%.5f', '%.5f'], delimiter="\t")
        np.savetxt("curve2.tsv", np.array([c[1] for c in curves]), fmt=['%.5f', '%.5f', '%.5f'], delimiter="\t")
        geotocreturn = sub.run([self.exe, "mode=dtw", "curve1=curve1.tsv", "curve2=curve2.tsv"], capture_output=True)
        twdensities = [float(line.split("=")[-1]) for line in geotocreturn.stdout.decode().split("\n")[3:]
                       if len(line) > 0 and "=" in line]
        os.chdir(originaldir)
        shutil.rmtree(pathlib.Path(tempdir))
        assert len(twdensities) == len(curves)
        return twdensities
    
    def analyze_topology_curves(self, curves: List) -> List[float]:
        originaldir = pathlib.Path.cwd()
        tempdir = tempfile.mkdtemp()
        os.chdir(tempdir)
        np.savetxt("curve1.tsv", np.array([c[0] for c in curves]), fmt=['%.5f', '%.5f', '%.5f'], delimiter="\t")
        np.savetxt("curve2.tsv", np.array([c[1] for c in curves]), fmt=['%.5f', '%.5f', '%.5f'], delimiter="\t")
        geotocreturn = sub.run([self.exe, "mode=all", "curve1=curve1.tsv", "curve2=curve2.tsv"], capture_output=True)
        tw, wr, lk, calccheck = [float(line.split("=")[-1]) for line in geotocreturn.stdout.decode().split("\n")[3:] 
                                 if len(line) > 0 and "=" in line]
        os.chdir(originaldir)
        shutil.rmtree(pathlib.Path(tempdir))
        # method returns: Twist, Writhe, Linking number
        return tw, wr, lk
    
    def analyze_topology_frames(self, bpframes: List) -> List[float]:
        originaldir = pathlib.Path.cwd()
        tempdir = tempfile.mkdtemp()
        os.chdir(tempdir)
        curve1 = np.array([np.array(frame[0]) for frame in bpframes])
        curve2 = np.array([np.array(frame[0])+np.array(frame[1][1])*10 for frame in bpframes])
        np.savetxt("curve1.tsv", curve1, fmt=['%.5f', '%.5f', '%.5f'], delimiter="\t")
        np.savetxt("curve2.tsv", curve2, fmt=['%.5f', '%.5f', '%.5f'], delimiter="\t")
        geotocreturn = sub.run([self.exe, "mode=all", "curve1=curve1.tsv", "curve2=curve2.tsv"], capture_output=True)
        tw, wr, lk, calccheck = [float(line.split("=")[-1]) for line in geotocreturn.stdout.decode().split("\n")[3:] 
                                 if len(line) > 0 and "=" in line]
        os.chdir(originaldir)
        shutil.rmtree(pathlib.Path(tempdir))
        # method returns: Twist, Writhe, Linking number
        return tw, wr, lk

