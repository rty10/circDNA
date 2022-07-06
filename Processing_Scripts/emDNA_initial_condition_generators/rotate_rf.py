#! /usr/bin/python
import os, math, shutil, subprocess, scipy
import numpy as np
PATH = os.getcwd()


def rf_reformat(filename, new_filename):
    name, ext = filename.split('.')
    infile =  open(filename, 'r')
    x = infile.readlines()
    infile.close()
    #os.system('rm '+filename)
    txt = [rf.rstrip('\n').split('#')[0] for rf in x]
    seq = []
    N_seq = int((txt[0].split()[0]))
    
    for i in range(0, len(txt)):
        if txt[i].startswith('...'):
            nt = txt[i].split()[2]
            if nt in basepairs:
                seq.append(nt)
            elif nt in basepair_fix:
                seq.append( bp_edit_dict[nt] )
    
    for i in range(0, len(seq)):
        j = 1+(i)*5
        refframe = np.array([[float(txt[j+2].split()[0]), float(txt[j+2].split()[1]), float(txt[j+2].split()[2])],
                             [float(txt[j+3].split()[0]), float(txt[j+3].split()[1]), float(txt[j+3].split()[2])],
                             [float(txt[j+4].split()[0]), float(txt[j+4].split()[1]), float(txt[j+4].split()[2])]], dtype='float64')
        txt[j] = '... ' +str(i+1) + ' ' +seq[i]+ ' ...'
        txt[j+1] = '{:>10} {:>10} {:>10}  # origin'.format(round(float(txt[j+1].split()[0]), 5), round(float(txt[j+1].split()[1]), 5), round(float(txt[j+1].split()[2]), 5) )
        txt[j+2] = '{:>10} {:>10} {:>10}  # x-axis'.format(round(refframe.item((0,0)), 5), round(refframe.item((0,1)), 5), round(refframe.item((0,2)), 5) )
        txt[j+3] = '{:>10} {:>10} {:>10}  # y-axis'.format(round(refframe.item((1,0)), 5), round(refframe.item((1,1)), 5), round(refframe.item((1,2)), 5) )
        txt[j+4] = '{:>10} {:>10} {:>10}  # z-axis'.format(round(refframe.item((2,0)), 5), round(refframe.item((2,1)), 5), round(refframe.item((2,2)), 5) )
    outfile = open(new_filename+'.'+ext, 'w')
    for line in txt:
        outfile.write(line+'\n')
    outfile.close()
    return


def Rotation_z(rot_angle):
    return np.array( [[np.cos(rot_angle), -1*np.sin(rot_angle), 0],
                      [np.sin(rot_angle), np.cos(rot_angle), 0],
                      [0, 0, 1]], dtype='float64')


def rf_rotate(filename, angle):
    """
    NOTE: Make sure to do this prior to generate circular reference frame file
    """
    infile =  open(filename, 'r')
    name = filename.split('.')[0]
    if '_std' in name:
        name = name.replace('_std', '')
    x = infile.readlines()
    infile.close()
    txt = [rf.strip('\n') for rf in x]
    seq = []
    N_seq = int((txt[0].split()[0]))
    
    for i in range(0, len(txt)):
        if txt[i].startswith('...'):
            seq.append(txt[i].split()[2])
    
    for x in range(1, 360//int(angle)):
        txt_rot = []
        txt_rot.append(txt[0])
        theta = np.pi*(float(x)/(1/float(2)*float(360//angle)))
        for i in range(0, len(seq)):
            rot = Rotation_z(theta)
            j = 1+(i)*5
            refframe = np.array([[float(txt[j+2].split()[0]), float(txt[j+2].split()[1]), float(txt[j+2].split()[2])],
                                 [float(txt[j+3].split()[0]), float(txt[j+3].split()[1]), float(txt[j+3].split()[2])],
                                 [float(txt[j+4].split()[0]), float(txt[j+4].split()[1]), float(txt[j+4].split()[2])]])
            
            refframe_rot = np.matmul(rot, refframe)            
            txt_rot.append(txt[j])
            txt_rot.append(txt[j+1])
            txt_rot.append('{:>10} {:>10} {:>10}  # x-axis'.format(round(refframe_rot.item((0,0)), 5), round(refframe_rot.item((0,1)), 5), round(refframe_rot.item((0,2)), 5) ) )
            txt_rot.append('{:>10} {:>10} {:>10}  # y-axis'.format(round(refframe_rot.item((1,0)), 5), round(refframe_rot.item((1,1)), 5), round(refframe_rot.item((1,2)), 5) ) )
            txt_rot.append('{:>10} {:>10} {:>10}  # z-axis'.format(round(refframe_rot.item((2,0)), 5), round(refframe_rot.item((2,1)), 5), round(refframe_rot.item((2,2)), 5) ) )

        outfile = open(name+'_r{:03d}.dat'.format(int(round(np.rad2deg(theta)))), 'w')
        for line in txt_rot:
            outfile.write(line+'\n')
        outfile.close()
        del txt_rot[:]
    return


def bpstep_insert(filename):
    """
    Function that will copy the first base pair coordinate frame in a reference frame document and paste it to the end.
    This is to be done on CIRCLES only.
    :param filename:
    :return:
    """
    f = open(filename, "r")
    name = filename.split('.')[0]
    ext = filename.split('.')[1]
    g = open(name+'_circ.'+ext, "w")
    x = f.readlines()
    f.close()
    
    N_seq = int((x[0].split()[0]))
    x[0] = '{:>5} base-pairs\n'.format(N_seq)
    i = 1
    for j in range(0, len(x)):
        g.write(x[j])
    while i < 6:
        g.write(x[i])
        i += 1
    g.close()
    return

