#! /usr/bin/python

import os

path = os.getcwd()
FILES = []
for filename in os.listdir(path):
    if filename.endswith('.dat'):
        FILES.append( filename.split('.')[0] )

for j in range(0, len(sorted(FILES))):
    FILENAME=FILES[j]
    infile = open(FILENAME+'.dat', 'r')
    data = infile.readlines()
    data = data[1:]
    size = int(len(data)/5)
    infile.close()    

    outfile = open("origins_"+FILENAME+".pdb", 'w')
    outfile.write("REMARK    PDB file from x3DNA reference frame origins from file: "+FILENAME+", October 2021, R Young\n")
    for i in range(0, size):
        bp = data[(5*i)].split()[2][0:1]
        X  = round(float(data[(5*i)+1].split()[0]), 3)
        Y  = round(float(data[(5*i)+1].split()[1]), 3)
        Z  = round(float(data[(5*i)+1].split()[2]), 3)
        outfile.write(("ATOM{:>7}  C1' {:>3} A{:>4}    {:>8}{:>8}{:>8}  1.00  1.00           C\n".format(str(i+1), "D"+bp, str(i+1), str(X), str(Y), str(Z))))
        del bp, X, Y, Z
    outfile.close()

    del data, size, FILENAME
del FILES

