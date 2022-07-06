#! /usr/bin/python
import os, sys
import math, scipy
import numpy as np
import pandas as pd
from scipy.spatial import distance
path = os.getcwd()
if "Young_Research" in path:
    initialpath = "C:\\Users\\Young_Research\\Documents\\Rutgers\\Research\\InitialConditions"
    ffpath = "C:\\Users\\Young_Research\\Documents\\Rutgers\\Research\\ForceFields"
else:
    initialpath = "/home/rty10/Documents/initialconditions"
    ffpath = "/home/rty10/Documents/forcefields"

ZERO = 0.00000001
ZERO2 = ZERO*ZERO
GLOBAL_O = np.array([0.0, 0.0, 0.0])
GLOBAL_X = np.array([1.0, 0.0, 0.0])
GLOBAL_Y = np.array([0.0, 1.0, 0.0])
GLOBAL_Z = np.array([0.0, 0.0, 1.0])
GLOBAL_FRAME = [GLOBAL_O, [GLOBAL_X, GLOBAL_Y, GLOBAL_Z]]


# ------- Functions to check for errors --------------------------
def _check_return_code(log_file):
    infile = open(log_file+".log", 'r')
    rfdata = infile.readlines()
    infile.close()
    for line in rfdata:
        if 'return code' in line:
            code = line.split()[2]
    return code


def _check_selfcontacts(filename, dataframe, starting_point, ending_point):
    '''
    filename = the reference frame file name 
    dataframe = needs to be a dataframe of (x, y, z) origin coordiates from the reference frame file
    starting_point = the index number to begin the collision check, in case of frozen steps (default = 0)
    ending_point   = the index number that ends region of collision check, in case of frozen steps (default = len(dataframe) )
    '''
    X = pd.DataFrame( distance.squareform( distance.pdist( dataframe ) ) )
    Y = X.where(np.triu(np.ones(X.shape)).astype(np.bool))
    Y = Y.stack().reset_index()
    Y.columns = ['bp_i','bp_j','distance']
    Y = Y[Y['distance']<=20.5].reset_index(drop=True)
    if starting_point != 0:
        Y = Y[(Y['bp_i']>=starting_point)&(Y['bp_j']>=starting_point)].reset_index(drop=True)
    if ending_point != len(X):
        Y = Y[(Y['bp_i']<=ending_point)&(Y['bp_j']<=ending_point)].reset_index(drop=True)
    contact_list = []
    if len(Y) > 0:
        for i in range(len(Y)):
            if abs(Y.at[i, 'bp_i'] - Y.at[i, 'bp_j']) > 11:
                contact_list.append((Y.at[i, 'bp_i'],Y.at[i, 'bp_j']))
        contact_list=sorted(list(set(contact_list)))
    del X, Y
    return contact_list
    

def _check_selfcontacts_constraints(filename, dataframe, constraintlength):
    X = pd.DataFrame( distance.squareform( distance.pdist( dataframe ) ) )
    contact_list = []
    for i in range(len(X)):
        for j in range(len(X)):
            if i > constraintlength and j > constraintlength:
                if X.at[i, j] <= 20.5 and abs(j-i) > 11:
                    contact_list.append(i)
    contact_list=sorted(list(set(contact_list)))
    del X
    return contact_list


def are_equal(f1, f2, threshhold=ZERO):
    return np.abs(np.array(f1) - np.array(f2)) < threshhold


# ------- Mathematics --------------------------
def norm_vec(vector):
    return vector/np.linalg.norm(vector)

def _vec_distance(vector):
    return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

def _plane_normal(point1, point2, point3):
    return norm_vec(np.cross(norm_vec(point2 - point1), norm_vec(point3 - point1)))
    
def _plane_projection(planenormal, vector):
    return norm_vec(vector - np.dot(vector, planenormal)*planenormal)

def _positions_angle(positionvector1, positionvector2, positionvector3):
    vec1 = norm_vec(positionvector1 - positionvector2)
    vec2 = norm_vec(positionvector3 - positionvector2)
    return np.rad2deg( np.arccos( np.dot(vec1, vec2) ) )

def _vecvec_angle(vector1, vector2):
    return np.rad2deg( np.arccos( np.dot( norm_vec(vector1), norm_vec(vector2) ) ) )

def _vecpln_angle(vector, planenormal):
    return np.rad2deg( np.arcsin( np.dot( norm_vec(vector), norm_vec(planenormal) ) ) )


def _loop_origins_omega_angles(firstconstraintnormal, lastconstraintnormal):
    vecA = -1*(norm_vec(firstconstraintnormal))
    vecB = norm_vec(lastconstraintnormal)
    return np.rad2deg( np.arccos( np.dot(vecA, vecB) ) )

def _loop_westcott_angles(dataframe_origins, dataframe_normals, loop_length, start_bp_id, midloop_bp_id, end_bp_id):
    '''
    '''
    start_origin = dataframe_origins.loc[start_bp_id].to_numpy()
    mid_origin   = dataframe_origins.loc[midloop_bp_id].to_numpy()
    end_origin   = dataframe_origins.loc[end_bp_id].to_numpy()
    
    start_norm = dataframe_normals.loc[start_bp_id].to_numpy()
    mid_norm   = dataframe_normals.loc[midloop_bp_id].to_numpy()
    end_norm   = dataframe_normals.loc[end_bp_id].to_numpy()
        
    enex_vector = start_origin - end_origin
    eta         = _vec_distance(enex_vector)/(3.4*loop_length)
    nhat        = _plane_normal(start_origin, mid_origin, end_origin)
    alpha       = _vecvec_angle(enex_vector, start_norm)
    beta        = _vecvec_angle(start_norm, _plane_projection(nhat, start_norm) )
    gamma       = _vecpln_angle(mid_norm, nhat)
    del start_origin, start_norm, mid_origin, mid_norm, end_origin, end_norm, enex_vector
    return eta, alpha, beta, gamma

def _loop_pinching(dataframe_origins, Ncirc, firstconstraint, lastconstraint):
    '''
    Get a square matrix of distances between the ith/jth origin pair along a loop.
    For each ith origin, determine the shortest jth distance.
    From the list of ith distances, save the minimal distance of this collection = "loop pinch point"
    '''
    minlist = {}
    INLOOP = 22
    INLOOPMAX = 3.4*INLOOP    
    df1 = dataframe_origins.copy()
    df2 = pd.DataFrame(distance.squareform( distance.pdist( df1[['x','y','z']] ) ), 
                       columns=np.arange(1, len(df1)+1), index=np.arange(1, len(df1)+1) )
    if int(lastconstraint) > int(firstconstraint):
        lst = [int(i) for i in np.arange(firstconstraint, lastconstraint+1)]
        df3 = df2.loc[[i for i in lst], [i for i in lst]]
        df3 = df3.reset_index(drop=False).rename(columns={'index':'bp'})
        df3.index+=1
    elif lastconstraint == 1:
        lst = [int(i) for i in np.arange(firstconstraint+1, Ncirc+1)]
        df3 = df2.loc[[i for i in lst], [i for i in lst]]
        df3 = df3.reset_index(drop=False).rename(columns={'index':'bp'})
        df3.index+=1
    else:
        lst = [int(i) for i in np.arange(1, lastconstraint)] + [int(i) for i in np.arange(firstconstraint+1, Ncirc+1)]
        df3 = df2.loc[[i for i in lst], [i for i in lst]]
        df3 = df3.reset_index(drop=False).rename(columns={'index':'bp'})
        df3a = df3.loc[:lastconstraint-2]
        df3b = df3.loc[lastconstraint-1:]
        df3  = pd.concat([df3b, df3a], ignore_index=True)
        df3.index+=1
        del df3a, df3b
    for i in lst:
        A    = [j for j in df3.index if df3.at[j, 'bp']==i]
        loop = [j for j in df3.index if A[0] <= j-INLOOP or A[0] >= j+INLOOP]
        X = df3.loc[[l for l in loop], i]
        minlist[i]=[df3.at[X.idxmin(), 'bp'], X.min()]
        del A, loop, X
    del df3

    A = np.array([minlist[i][1] for i in minlist])
    for key, value in minlist.items():
        if A.min() == value[1] and key < value[0]:
            pinchpoints = "(" + str(key) + ", " + str(value[0]) + ")"
            pinchdistance = value[1]
            pinchbpdis    = pinchdistance/3.4
    del A, minlist, INLOOP, INLOOPMAX, lst
    del df1, df2
    return pinchpoints, pinchdistance, pinchbpdis


def _ellipticity(pca_dataframe):
    return (pca_dataframe.e1.max() - pca_dataframe.e1.min())/(pca_dataframe.e2.max() - pca_dataframe.e2.min())


# ------- DATAFRAME FUNCTIONS --------------------------
def refframe_get_sequence(refframe_filename):
    '''
    '''
    BPDICT = {"A-T":"A", "C-G":"C", "G-C":"G", "T-A":"T"}
    SEQ = ''
    infile = open(refframe_filename, 'r')
    rfdata = infile.readlines()
    infile.close()
    rfdata = [i.rstrip('\n').replace('  # origin','').replace('  # x-axis','').replace('  # y-axis','').replace('  # z-axis','').split() for i in rfdata]
    N = int(rfdata[0][0])
    for j in range(0, N):
        BASE = BPDICT[rfdata[5*j+1][2]]
        SEQ += BASE
        del BASE
    del BPDICT
    return SEQ


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

def df_read_refframe_normals(refframe_filename):
    '''
    '''
    df   = pd.DataFrame(columns=['x','y','z'])
    infile = open(refframe_filename, 'r')
    rfdata = infile.readlines()
    infile.close()
    rfdata = [i.rstrip('\n').replace('  # origin','').replace('  # x-axis','').replace('  # y-axis','').replace('  # z-axis','').split() for i in rfdata]
    N = int(rfdata[0][0])
    for j in range(0, N):
        df.at[j, ['x','y','z']]  = rfdata[ 5*j + 5]
    df = df.astype('float')
    df.index += 1
    return df

def _gyration_origindf(dataframe):
    df = dataframe.copy()
    df['rg'] = (dataframe.x-dataframe.x.mean())**2 + (dataframe.y-dataframe.y.mean())**2 + (dataframe.z-dataframe.z.mean())**2
    X = np.sqrt( (df.rg).mean() )
    del df
    return X

def _pca_origindf(refframe_dataframe):
    df = refframe_dataframe.copy()
    N = len(df)
    df['vx'] = df.x - df.x.mean()
    df['vy'] = df.y - df.y.mean()
    df['vz'] = df.z - df.z.mean()

    cxx = ( (df.vx*df.vx).sum() - (df.vx.sum() * df.vx.sum()) ) / (N-1)
    cxy = ( (df.vx*df.vy).sum() - (df.vx.sum() * df.vy.sum()) ) / (N-1)
    cxz = ( (df.vx*df.vz).sum() - (df.vx.sum() * df.vz.sum()) ) / (N-1)
    cyx = ( (df.vy*df.vx).sum() - (df.vy.sum() * df.vx.sum()) ) / (N-1)
    cyy = ( (df.vy*df.vy).sum() - (df.vy.sum() * df.vy.sum()) ) / (N-1)
    cyz = ( (df.vy*df.vz).sum() - (df.vy.sum() * df.vz.sum()) ) / (N-1)
    czx = ( (df.vz*df.vx).sum() - (df.vz.sum() * df.vx.sum()) ) / (N-1)
    czy = ( (df.vz*df.vy).sum() - (df.vz.sum() * df.vy.sum()) ) / (N-1)
    czz = ( (df.vz*df.vz).sum() - (df.vz.sum() * df.vz.sum()) ) / (N-1)

    covar = np.matrix([ [cxx, cxy, cxz], [cyx, cyy, cyz], [czx, czy, czz] ])
    evals, evecs = np.linalg.eig( covar )
    idx = np.argsort(evals)[::-1][: len(evals)]
    evals = evals[idx]
    evecs = evecs[:, idx]

    adj = pd.DataFrame(np.dot( df[['vx', 'vy', 'vz']], evecs), columns=['e1','e2','e3'])
    del df, covar, evals, evecs, idx
    return adj
    

def df_read_bpsteppars(par_filename, bend=False):
    '''
    This function loads an optimized parameter file, calculates the energy per base-pair step and the dimer bend angles
    Returns a parameter dataframe
    '''
    df = pd.read_csv(par_filename, sep='\s+', skiprows=2)
    df.columns = df.columns.str.lower()
    df = df.drop(columns=['shear','stretch','stagger','buckle','prop-tw','opening'])
    for i in df.index:
        if i > 0:
            df.at[i, 'dimer'] = df.at[i-1, '#'][0:1] + df.at[i, '#'][0:1]            
    if bend==True:
        df['bend'] = np.sqrt(df['tilt']**2 + df['roll']**2)
        df = df[['#','dimer','tilt','roll','bend','twist','shift','slide','rise']]
    else:
        df = df[['#','dimer','tilt','roll','twist','shift','slide','rise']]
    return df

def _par_tetrameric(bps_dataframe):
    '''
    This function adds a column of tetramer sequence steps to a base-pair step parameter dataframe
    '''
    df = bps_dataframe.copy()
    SEQ=''.join([df.at[i,'#'][0:1] for i in range(len(df))])
    SEQ = 'x'+SEQ+'x'
    for i in df.index:
        if i > 0:
            df.at[i, 'tetramer'] = SEQ[i-1:i+3]
    del SEQ    
    df = df[['#','dimer','tetramer','tilt','roll','twist','shift','slide','rise']]
    return df

def _par_reststate_analysis(sequence, dfrs):
    '''
    This function generates a parameter dataframe of a sequence-specific rest state.
    This function requires a sequence.
    Returns a parameter dataframe
    '''
    
    df = pd.DataFrame(columns=['#','Tilt','Roll','Twist','Shift','Slide','Rise'], index=[i for i in range(len(sequence))])
    bases = {'A':'A-T','T':'T-A','C':'C-G','G':'G-C'}
    for i in df.index:
        df.at[i, '#'] = bases[sequence[i]]
    del bases
    if 'step' not in df.columns.to_list():
        for i in df.index:
            if i > 0:
                df.at[i, 'dimer'] = df.at[i-1, '#'][0:1] + df.at[i, '#'][0:1]
    for i in df.index:
        STEP = df.loc[i]['dimer']
        for PAR in ['Tilt','Roll','Twist','Shift','Slide','Rise']:
            if i == 0 :
                df.at[i, PAR]=0
            else:
                df.at[i, PAR]=dfrs.loc[STEP][PAR]
        del STEP
    return df    


def _par_energetic_analysis(dataframe, dfrs, dffc, tet=False):
    '''
    This function loads an optimized parameter file, calculates the energy per base-pair step and the dimer bend angles
    Returns a parameter dataframe
    '''
    df = dataframe.copy()
    if 'dimer' not in df.columns.to_list():
        for i in df.index:
            df.at[i, 'dimer'] = df.at[i-1, '#'][0:1] + df.at[i, '#'][0:1]
            
    if tet==True:
        df['step']=df['tetramer']
    else:
        df['step']=df['dimer']
    
    for i in df.index:
        if i > 0:
            vector = df.loc[i]
            parvec = dfrs.loc[vector['step']].to_numpy()
            fcmat  = dffc.loc[vector['step']].to_numpy().reshape((6,6))
            dimvec = vector[['tilt','roll','twist','shift','slide','rise']].to_numpy()
            diffvec = dimvec-parvec
            energy = 1/2*np.transpose(diffvec).dot( fcmat.dot(diffvec) )
            df.at[i, 'energy']=energy
            del vector, parvec, fcmat, dimvec, diffvec, energy
    df = df.drop('step', axis=1)
    return df    

def df_read_topology(topology_filename):
    '''
    '''
    infile = open(topology_filename, 'r')
    topo_file_list = infile.readlines()
    topo_file_list = [i.rstrip('\n') for i in topo_file_list]
    infile.close()
    topo_file_list = topo_file_list[-4:]
    for i in range(0, len(topo_file_list)):
        if 'Wr ' in topo_file_list[i]:
            wr = float( topo_file_list[i].split('=')[1] )
        elif 'Tw ' in topo_file_list[i]:
            tw = float( topo_file_list[i].split('=')[1] )
        elif 'Lk ' in topo_file_list[i]:
            lk = int( topo_file_list[i].split('=')[1] )
    del topo_file_list[:]
    return wr, tw, lk
    
def df_read_pdb(pdb_filename):
    '''
    '''
    infile = open(pdb_filename, 'r')
    topo_file_list = infile.readlines()
    infile.close()
    topo_file_list = [i.rstrip('\n') for i in topo_file_list if 'ATOM' in i]
    pdb_header=['record','atom_serial','atom_name','loc_ind','res_name','chain_id','res_seq_num','insert_code','x','y','z','occ','beta','segment_id','element']
    df=pd.DataFrame(columns=pdb_header, index=[i for i in range(len(topo_file_list))])
    for i in range(len(topo_file_list)):
        df.at[i, 'record']=topo_file_list[i][0:4]
        df.at[i, 'atom_serial']=topo_file_list[i][6:11]
        df.at[i, 'atom_name']=topo_file_list[i][12:16]
        df.at[i, 'loc_ind']=topo_file_list[i][16:17]
        df.at[i, 'res_name']=topo_file_list[i][17:20]
        df.at[i, 'chain_id']=topo_file_list[i][21:22]
        df.at[i, 'res_seq_num']=topo_file_list[i][22:26]
        df.at[i, 'insert_code']=topo_file_list[i][26:27]
        df.at[i, 'x']=topo_file_list[i][30:38]
        df.at[i, 'y']=topo_file_list[i][38:46]
        df.at[i, 'z']=topo_file_list[i][46:54]
        df.at[i, 'occ']=topo_file_list[i][54:60]
        df.at[i, 'beta']=topo_file_list[i][60:66]
        df.at[i, 'segment_id']=topo_file_list[i][72:76]
        df.at[i, 'element']=topo_file_list[i][76:78]
    df[['x','y','z']]=df[['x','y','z']].astype(float)
    df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    return df


# --!-- Dataframe: loading forcefield values ---------------------------------------------------------------
def df_read_intrinsic_state(forcefieldpath, forcefield):

    if os.path.isfile(forcefieldpath+'/RestStateParameters/StepParameters_'+forcefield+'.txt') == True:
        infile = open(forcefieldpath+'/RestStateParameters/StepParameters_'+forcefield+'.txt', 'r')
    elif os.path.isfile(forcefieldpath+'/RestStateParameters/StepParameters_'+forcefield.capitalize()+'.txt') == True:
        infile = open(forcefieldpath+'/RestStateParameters/StepParameters_'+forcefield.capitalize()+'.txt', 'r')    

    indata = infile.readlines()
    indata = [i.replace("={"," ").replace(", "," ").replace("}","").rstrip('\n').split() for i in indata]
    header = ['step','Tilt','Roll','Twist','Shift','Slide','Rise']
    for i in range(0, len(indata)):
        for j, x in enumerate(indata[i]):
            try:
                indata[i][j] = float(x)
            except ValueError:
                pass
    df1 = pd.DataFrame.from_records(indata, columns=header)
    df1 = df1.set_index('step')
    del indata
    return df1

# --- load force constants 
def df_read_force_constants(forcefieldpath, forcefield):
    
    if os.path.isfile(forcefieldpath+'/ForceConstants/ForceConstants_'+forcefield+'.txt') == True:
        infile = open(forcefieldpath+'/ForceConstants/ForceConstants_'+forcefield+'.txt', 'r')
    elif os.path.isfile(forcefieldpath+'/ForceConstants/ForceConstants_'+forcefield.capitalize()+'.txt') == True:
        infile = open(forcefieldpath+'/ForceConstants/ForceConstants_'+forcefield.capitalize()+'.txt', 'r')    
    else:
        infile = open(forcefieldpath+'/ForceConstants/ForceConstants_IdealDNA.txt', 'r')
        
    indata2 = infile.readlines()
    indata2 = [i.replace("={"," ").replace(", "," ").replace("}","").rstrip('\n').split() for i in indata2]
    for i in range(0, len(indata2)):
        for j, x in enumerate(indata2[i]):
            try:    indata2[i][j] = float(x)
            except ValueError:    pass
    header2 =['step',
    'TiltTilt','TiltRoll','TiltTwist','TiltShift','TiltSlide','TiltRise',
    'RollTilt','RollRoll','RollTwist','RollShift','RollSlide','RollRise',
    'TwistTilt','TwistRoll','TwistTwist','TwistShift','TwistSlide','TwistRise',
    'ShiftTilt','ShiftRoll','ShiftTwist','ShiftShift','ShiftSlide','ShiftRise',
    'SlideTilt','SlideRoll','SlideTwist','SlideShift','SlideSlide','SlideRise',
    'RiseTilt','RiseRoll','RiseTwist','RiseShift','RiseSlide','RiseRise']
    df2 = pd.DataFrame.from_records(indata2, columns=header2)
    df2 = df2.set_index('step')
    del indata2
    return df2


