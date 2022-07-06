#! /usr/bin/python

import os, math, Bio
import pandas as pd
import numpy as np
from Bio.PDB import *

ZERO = 0.000000000000001
ZERO2 = ZERO*ZERO
GLOBAL_O = np.array([0.0, 0.0, 0.0])
GLOBAL_X = np.array([1.0, 0.0, 0.0])
GLOBAL_Y = np.array([0.0, 1.0, 0.0])
GLOBAL_Z = np.array([0.0, 0.0, 1.0])
GLOBAL_FRAME = [GLOBAL_O, [GLOBAL_X, GLOBAL_Y, GLOBAL_Z]]

def np_read_refframe(refframe_filename, normalize=False):
    '''
    Load refframe file. Convert each 5-line frame into a single line, 2-dim array
    '''
    ref_frames = []
    infile = open(refframe_filename, "r")
    data = infile.readlines()
    for _i in range(int(data[0].split()[0])):
        description = data[_i * 5 + 1].strip()
        origin = [float(f) for f in data[_i * 5 + 2].split()[:3]]
        x = [float(f) for f in data[_i * 5 + 3].split()[:3]]
        y = [float(f) for f in data[_i * 5 + 4].split()[:3]]
        z = [float(f) for f in data[_i * 5 + 5].split()[:3]]
        if normalize:
            x = (x / np.linalg.norm(x)).tolist()
            y = (y / np.linalg.norm(y)).tolist()
            z = (z / np.linalg.norm(z)).tolist()
        ref_frames.append([description, origin, [x, y, z]])
    infile.close()
    return ref_frames

def _norm_vector(vector):
    return vector / np.sqrt( np.dot(vector, vector) )

def _basis_frame(frame):
    '''
    The reference frame file is formatted such that the origin and bases are a 4x3 matrix.
    Basis vectors are meant to be columns not rowse
    Proper form: B_i = [[O], [d_1, d_2, d_3]]_i
    '''
    return np.transpose(frame)

def _axis_angle_rotation(angle_radians, unit_vector):
    '''
    General rotation matrix, R(phi, u).
    To rotate vector p about some unit vector u through angle phi to new position q.
    Return the rotation matrix that will be used with vector p
    '''
    A = angle_radians
    U = _norm_vector(unit_vector)
    matrix_1 = np.array([np.cos(A)+(U[0]*U[0])*(1-np.cos(A)),      (U[0]*U[1])*(1-np.cos(A))-U[2]*np.sin(A), (U[0]*U[2])*(1-np.cos(A))+U[1]*np.sin(A)])
    matrix_2 = np.array([(U[0]*U[1])*(1-np.cos(A))+U[2]*np.sin(A), np.cos(A)+(U[1]*U[1])*(1-np.cos(A)),      (U[1]*U[2])*(1-np.cos(A))-U[0]*np.sin(A)])
    matrix_3 = np.array([(U[0]*U[2])*(1-np.cos(A))-U[1]*np.sin(A), (U[1]*U[2])*(1-np.cos(A))+U[0]*np.sin(A), np.cos(A)+(U[2]*U[2])*(1-np.cos(A))])
    ROT = np.array([[matrix_1], [matrix_2], [matrix_3]])
    del matrix_1, matrix_2, matrix_3
    return ROT
    
def _rotation_matrix_byframes(mobile_frame, target_frame):
    return np.matmul(np.transpose(mobile_frame), target_frame)

def _rotation_matrix_bybases(mobile_basis_frame, target_basis_frame):
    return np.matmul(np.linalg.inv(target_basis_frame), mobile_basis_frame)

def _transformation_by_frames_rotation(coordinate_array, rotation_matrix):
    return np.dot(coordinate_array, rotation_matrix)

def _transformation_by_bases_rotation(coordinate_array, rotation_matrix):
    return np.dot(coordinate_array, np.transpose(rotation_matrix))


def constructed_basis_frame(coord_array, coord_ref_pt):   
    reference_vector_00 = coord_array[coord_ref_pt]
    reference_vector_b1 = coord_array[coord_ref_pt - 1]
    reference_vector_f1 = coord_array[coord_ref_pt + 1]
    
    vec_1  = reference_vector_f1 - reference_vector_00
    nvec_1 = _norm_vector( vec_1 )
    tvec = reference_vector_00 - reference_vector_b1
    ntvec = _norm_vector( tvec )
    vec_3 = np.cross(nvec_1, ntvec)
    nvec_3 = _norm_vector( vec_3 )
    vec_2 = np.cross(nvec_3, nvec_1)
    nvec_2 = _norm_vector( vec_2 )
    
    MATRIX = np.array( [nvec_1, nvec_2, nvec_3] )
    
    del reference_vector_00, reference_vector_b1, reference_vector_f1
    del vec_1, tvec, vec_2, vec_3, nvec_1, ntvec, nvec_2, nvec_3
    return MATRIX
