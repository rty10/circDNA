
import os, sys, shutil, subprocess, matplotlib
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import spatial
path = os.getcwd()

# Defined Functions ---------------------------------------------------------------------------------------------------------------

def gyration(filename):
    df     = pd.DataFrame(columns=['x','y','z'])
    infile = open(filename, 'r')
    rfdata = infile.readlines()
    infile.close()
    rfdata = [i.rstrip('\n').split() for i in rfdata]
    N = int(rfdata[0][0])
    for j in range(0, N):
        df.at[j, ['x','y','z']]  = rfdata[ 5*j + 2]
    df = df.astype('float')
    df['rg'] = (df.x-df.x.mean())**2 + (df.y-df.y.mean())**2 + (df.z-df.z.mean())**2
    R = np.sqrt( (df.rg).mean() )
    del rfdata, df
    return R


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
                               "Dx","Dy","Dz",
                               "-d1x","-d1y","-d1z",
                               "Dcenter","theta",
                               "px","py","pz","qx","qy","qz"])    

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
        orn = np.array(df.loc[i, ["Ox","Oy","Oz"]])
        ocr = np.array(df.loc[i, ["-d1x","-d1y","-d1z"]])
        
        df.loc[i, ["Dx","Dy","Dz"]] = com - orn
        df.loc[i, "Dcenter"]        = sp.spatial.distance.euclidean(com, orn)
        dc = (com - orn)/sp.spatial.distance.euclidean(com, orn)
        
        df.loc[i, "theta"] = np.degrees( np.arccos( np.dot(dc, ocr) ) )
        
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
    elif len(Pdata) == 2*Ncirc:
        for i in range(0, len(Odata)):
            df.loc[i, ["px","py","pz"]] = float(Pdata[i][0]), float(Pdata[i][1]), float(Pdata[i][2])
            j = ((2*Ncirc)-i)-1
            df.loc[i, ["qx","qy","qz"]] = float(Pdata[j][0]), float(Pdata[j][1]), float(Pdata[j][2])
    del rfdata, PDBdata, Odata, Xdata, Pdata
    #df = circ_minor_groove(df)
    #df = circ_major_groove(df)
    return df, com


def load_opt_par_dataframe(filepath, Nseq):
    infile = open(filepath, 'r')
    indata = infile.readlines()
    infile.close()
    indata = [i.rstrip('\n').split() for i in indata]
    header = ['basepair',
    'Shear','Stretch','Stagger','Buckle','Prop-Tw','Opening',
    'Shift','Slide','Rise','Tilt','Roll','Twist']
    indata = indata[3:]
    for i in range(0, len(indata)):
        for j, x in enumerate(indata[i]):
            try:
                indata[i][j] = float(x)
            except ValueError:
                pass
    df = pd.DataFrame.from_records(indata, columns=header)
    del indata
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

def load_forcefield(ffpath, forcefield):
    """
    Function to generate a dataframe with rest state values from an optimization forcefield
    and another dataframe with elastic force constants of the same forcefield
    """
    if 'rt' not in forcefield:
        for filename in os.listdir(ffpath+'/RestStateParameters'):
            if forcefield in filename:
                infile1 = open(ffpath+'/RestStateParameters/'+filename, 'r')
        if forcefield == 'Olson1998' or forcefield == 'olson':
            infile2 = open(ffpath+'/ForceConstants/ForceConstants_Olson1998.txt', 'r')
        else:
            infile2 = open(ffpath+'/ForceConstants/ForceConstants_IdealDNA.txt', 'r')
    else:
        ff_name = forcefield.split('-')[0]
        for filename in os.listdir(ffpath+'/RestStateParameters'):
            if ff_name in filename:
                infile1 = open(ffpath+'/RestStateParameters/'+filename, 'r')
        infile2 = open(ffpath+'/ForceConstants/ForceConstants_Coleman2003.txt', 'r')
    
            
    indata1 = infile1.readlines()
    indata1 = [i.replace("={"," ").replace(", "," ").replace("}","").rstrip('\n').split() for i in indata1]
    for i in range(0, len(indata1)):
        for j, x in enumerate(indata1[i]):
            try:	indata1[i][j] = float(x)
            except ValueError:	pass
    indata2 = infile2.readlines()
    indata2 = [i.replace("={"," ").replace(", "," ").replace("}","").rstrip('\n').split() for i in indata2]
    for i in range(0, len(indata2)):
        for j, x in enumerate(indata2[i]):
            try:	indata2[i][j] = float(x)
            except ValueError:	pass
    header1 = ['dimer','tilt','roll','twist','shift','slide','rise']
    header2 =['dimer',
    'TiltTilt','TiltRoll','TiltTwist','TiltShift','TiltSlide','TiltRise',
    'RollTilt','RollRoll','RollTwist','RollShift','RollSlide','RollRise',
    'TwistTilt','TwistRoll','TwistTwist','TwistShift','TwistSlide','TwistRise',
    'ShiftTilt','ShiftRoll','ShiftTwist','ShiftShift','ShiftSlide','ShiftRise',
    'SlideTilt','SlideRoll','SlideTwist','SlideShift','SlideSlide','SlideRise',
    'RiseTilt','RiseRoll','RiseTwist','RiseShift','RiseSlide','RiseRise']
    df1 = pd.DataFrame.from_records(indata1, columns=header1)
    df2 = pd.DataFrame.from_records(indata2, columns=header2)
    indata2 = infile2.readlines()
    del indata1, indata2
    return df1, df2

def insert_bpstep_seq_circular(opt_par_dataframe):
    """
    From a dataframe with a column of base-pairs, generate the dimer and tetramer for each bp.
    Add new dimer and tetramer columns to dataframe.
    Note: this is for circular constructions
    """
    bpseq  = opt_par_dataframe['basepair']
    for k in range(0, len(bpseq)):
        first, second, third, fourth = k-2, k-1, k, k+1 
        if k-2 == -2:
            first  = (Nseq - 2)
            second = (Nseq - 1)
        elif k-2 == -1:
            first = (Nseq - 1)
        elif k == Nseq-1:
            fourth = 0
        elif k == Nseq:
            second = (Nseq-1)
            third  = 0
            fourth = 1
        a = bpseq[first].split('-')[0]
        b = bpseq[second].split('-')[0]
        c = bpseq[third].split('-')[0]
        d = bpseq[fourth].split('-')[0]
        dimerstep = "".join((b, c))
        tetrastep = "".join((a, b, c, d))
        opt_par_dataframe.at[k, 'dimer'] = dimerstep
        opt_par_dataframe.at[k, 'tetramer'] = tetrastep
    return opt_par_dataframe

def insert_bps_bend(opt_par_dataframe):
    """
    Function that takes the tilt and roll columns from a loaded dataframe and determines the bend angle.
    """
    for k in range(0, len(opt_par_dataframe)):
        x = float(opt_par_dataframe.loc[k, 'Tilt'])
        y = float(opt_par_dataframe.loc[k, 'Roll'])
        opt_par_dataframe.loc[k, 'Bend'] = float(np.sqrt(x**2 + y**2))
    return opt_par_dataframe

def insert_bps_energy(Nseq, opt_par_dataframe, reststate_par_dataframe, elastic_constants_dataframe):
    """
    Function that determines the energy per base-pair step.
    Must have:
    - loaded rest state dataframe
    - loaded elastic force constant dataframe
    - column with dimer and/or tetramer steps
    """
    opt_par_dataframe.loc[0, 'Energy'] = float(0)
    for k in range(1, len(opt_par_dataframe)):
        dim = opt_par_dataframe.loc[k, 'dimer']
        oshift, oslide, orise, otilt, oroll, otwist = [j for j in opt_par_dataframe.loc[k, 'Shift':'Twist']]
        A = np.array([otilt, oroll, otwist, oshift, oslide, orise])
        for j in range(0, len(reststate_par_dataframe)):
            if reststate_par_dataframe.loc[j, 'dimer'] == dim:
                B = np.array([z for z in reststate_par_dataframe.loc[j, 'tilt':'rise']])
        for j in range(0, len(elastic_constants_dataframe)):
            if elastic_constants_dataframe.loc[j, 'dimer'] == dim:
                F = np.array([[z for z in elastic_constants_dataframe.loc[j, 'TiltTilt':'TiltRise']],
                              [z for z in elastic_constants_dataframe.loc[j, 'RollTilt':'RollRise']],
                              [z for z in elastic_constants_dataframe.loc[j, 'TwistTilt':'TwistRise']],
                              [z for z in elastic_constants_dataframe.loc[j, 'ShiftTilt':'ShiftRise']],
                              [z for z in elastic_constants_dataframe.loc[j, 'SlideTilt':'SlideRise']],
                              [z for z in elastic_constants_dataframe.loc[j, 'RiseTilt':'RiseRise']]])
        opt_par_dataframe.loc[k, 'Energy'] = (1/2)*np.trace( (A-B) * F * (A-B) )
    return opt_par_dataframe


# --- Output new parameter file ---
def newfile_bpsdata(Nseq, main_dataframe, outputfilepath, outputname):
    """
    bp, dimer, tetramer, tilt, roll, bend, twist, energy, dcenter, anglecenter, Wmaj, Wmin, ...
    """
    A = main_dataframe
    outfile = open(outputfilepath+'/'+outputname+'_bps-data.txt', 'w')
    outfile.write(str(Nseq)+'  # base pairs\n')
    if not "initial" in str(outputfilepath):
        outfile.write("{:<4}{:>6}{:>9}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(
        'bp','dimer','tetramer','Tilt','Roll','Bend','Twist','Energy','Dcenter','theta','W-maj','W-min'
        ))
        for i in range(0, len(A)):
            outfile.write("{:<4}{:>6}{:>9}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(
            A.loc[i, 'basepair'], A.loc[i, 'dimer'], A.loc[i, 'tetramer'],
            round(A.loc[i, 'Tilt'],5),round(A.loc[i, 'Roll'],5),round(A.loc[i, 'Bend'],5),round(A.loc[i, 'Twist'],5),round(A.loc[i, 'Energy'], 5),
            round(A.loc[i, 'Dcenter'], 5), round(A.loc[i, 'theta'], 5), round(A.loc[i, 'W-maj'], 5), round(A.loc[i, 'W-min'], 5) 
            ))
        outfile.close()
    else:
        outfile.write("{:<4}{:>6}{:>9}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(
        'bp','dimer','tetramer','Tilt','Roll','Bend','Twist','Dcenter','theta','W-maj','W-min'
        ))
        for i in range(0, len(A)):
            outfile.write("{:<4}{:>6}{:>9}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(
            A.loc[i, 'basepair'], A.loc[i, 'dimer'], A.loc[i, 'tetramer'],
            round(A.loc[i, 'Tilt'],5),round(A.loc[i, 'Roll'],5),round(A.loc[i, 'Bend'],5),round(A.loc[i, 'Twist'],5),
            round(A.loc[i, 'Dcenter'], 5), round(A.loc[i, 'theta'], 5), round(A.loc[i, 'W-maj'], 5), round(A.loc[i, 'W-min'], 5) 
            ))
        outfile.close()
    return

def load_optdetailed_df(filepath):
    infile = open(filepath, 'r')
    indata = infile.readlines()
    infile.close()
    indata = [i.rstrip('\n').split() for i in indata]
    for i in range(0, len(indata)):
        for j, x in enumerate(indata[i]):
            try:
                indata[i][j] = float(x)
            except ValueError:
                pass
    df = pd.DataFrame.from_records(indata[2:], columns=indata[1:2])
    df.index = [i for i in range(1, len(df)+1)]
    del indata
    return df

def df_reststate(path):
    infile1 = open(path, 'r')
    indata1 = infile1.readlines()
    infile1.close()            
    indata1 = [i.replace("={"," ").replace(", "," ").replace("}","").rstrip('\n').split() for i in indata1]
    rsdf = pd.DataFrame.from_records(indata1, columns=['dimer','Tilt','Roll','Twist','shift','slide','rise'])
    rsdf = rsdf.astype({'Tilt':"float64",'Roll':"float64",'Twist':"float64",'shift':"float64",'slide':"float64",'rise':"float64"})
    rsdf = rsdf.set_index('dimer')
    return rsdf
