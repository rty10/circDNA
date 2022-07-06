#! /usr/bin/python
import os, shutil, subprocess, scipy
import sys
import numpy as np
import pandas as pd
# --- Load force field parameters 
def load_forcefield(ffpath, forcefield):
	for filename in os.listdir(ffpath+'/RestStateParameters'):
		if forcefield in filename:
			infile1 = open(ffpath+'/RestStateParameters/'+filename, 'r')
	indata1 = infile1.readlines()
	indata1 = [i.replace("={"," ").replace(", "," ").replace("}","").rstrip('\n').split() for i in indata1]
	for i in range(0, len(indata1)):
		for j, x in enumerate(indata1[i]):
			try:	indata1[i][j] = float(x)
			except ValueError:	pass
	if forcefield == 'Olson1998' or forcefield == 'olson':
		infile2 = open(ffpath+'/ForceConstants/ForceConstants_Olson1998.txt', 'r')
	else:
		infile2 = open(ffpath+'/ForceConstants/ForceConstants_IdealDNA.txt', 'r')
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

# --- New data for optimized parameter data

def insert_bpstep_seq(Nseq, opt_par_dataframe):
	bpseq  = opt_par_dataframe['basepair']
	for k in range(0, (Nseq+1)):
		first, second, third, fourth = k-2, k-1, k, k+1 
		if k-2 == -2:
			first  = (Nseq - 2)
			second = (Nseq - 1)
		elif k-2 == -1:
			first = (Nseq - 1)
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
	for k in range(0, len(opt_par_dataframe)):
		x = float(opt_par_dataframe.loc[k, 'Tilt'])
		y = float(opt_par_dataframe.loc[k, 'Roll'])
		opt_par_dataframe.loc[k, 'Bend'] = float(np.sqrt(x**2 + y**2))
	return opt_par_dataframe

def insert_bps_energy(Nseq, opt_par_dataframe, reststate_par_dataframe, elastic_constants_dataframe):
	opt_par_dataframe.loc[0, 'energy'] = float(0)
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
		opt_par_dataframe.loc[k, 'energy'] = (1/2)*np.trace( (A-B) * F * (A-B) )
	return opt_par_dataframe

# --- Load optimized parameter file and insert new data
# -!- WARNING: Calling functions within functions
def load_opt_par_dataframe(filepath, Nseq, reststate_par_dataframe, elastic_constants_dataframe):
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
	df = insert_bps_bend(df)
	df = insert_bpstep_seq(Nseq, df)
	df = insert_bps_energy(Nseq, df, reststate_par_dataframe, elastic_constants_dataframe)
	return df


# --- Output new parameter file
def new_par_file(Nseq, opt_par_dataframe, outputfilepath, outputname):
	A = opt_par_dataframe
	outfile = open(outputfilepath+'/'+outputname+'.par', 'w')
	outfile.write(str(Nseq)+'  # base pairs\n')
	outfile.write('0 # ***local base-pair & step parameters***\n')
	outfile.write("{:<6} {:<6} {:<6} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}  {:>10} {:>10}\n".format(
	'#basepair', 'dimer','tetramer','Shear','Stretch','Stagger','Buckle','Prop-Tw','Opening',
	'Shift','Slide','Rise','Tilt','Roll','Twist','Bend','Energy'
	))
	for i in range(0, len(A)):
		outfile.write("{:<6} {:<6} {:<6} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}  {:>10} {:>10}\n".format(
		A.loc[i, 'basepair'], A.loc[i, 'dimer'], A.loc[i, 'tetramer'],
		round(A.loc[i, 'Shear'],5), round(A.loc[i, 'Stretch'],5), round(A.loc[i, 'Stagger'],5), round(A.loc[i, 'Buckle'],5), round(A.loc[i, 'Prop-Tw'],5), round(A.loc[i, 'Opening'],5),
		round(A.loc[i, 'Shift'],5), round(A.loc[i, 'Slide'],5),   round(A.loc[i, 'Rise'],5),    round(A.loc[i, 'Tilt'],5),   round(A.loc[i, 'Roll'],5),    round(A.loc[i, 'Twist'],5),
		round(A.loc[i, 'Bend'],5),  round(A.loc[i, 'energy'], 5)
		))
	outfile.close()
	return


