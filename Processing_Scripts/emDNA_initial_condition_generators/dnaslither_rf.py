#! /usr/bin/python
import os, math, shutil, subprocess, scipy
import numpy as np
PATH = os.getcwd()


def rf_reformat(filename, new_filename):
    name = filename.split('.')[0]
    ext = filename.split('.')[1]
    infile =  open(filename, 'r')
    x = infile.readlines()
    infile.close()
    #os.system('rm '+filename)
    txt = [rf.strip('\n').split('#')[0] for rf in x]
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

def rf_slither(filename, slither_count):
    """
    NOTE: Make sure to do this prior to generate circular reference frame file
    """
    infile =  open(filename, 'r')
    name = filename.split('.')[0]
    x = infile.readlines()
    infile.close()
    if '_std' in name:
        name = name.replace('_std', '')
    txt = [rf.strip('\n') for rf in x]
    seq = []
    N_seq = int((txt[0].lstrip().split(' ')[0]))
    for i in range(1, N_seq+1):
        j = 1+(i-1)*5
        seq.append(txt[j].split(' ')[2])
    slseq = seq
    for i in range(1, slither_count):
        slseq = slseq[-1:] + slseq[:-1]
        for j in range(1, N_seq+2):
            k = 1+(j-1)*5
            txt[k] = '... ' +str(j) + ' ' +slseq[j-1]+ ' ...'
        outfile = open(name+'_s{:03d}'.format(i)+'.dat', 'w')
        for line in txt:
            outfile.write(line+'\n')
        outfile.close()
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

