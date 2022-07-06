#! /usr/bin/python
import argparse
'''
Robert Young, Rutgers University, 2021

Python3 script that loads an emDNA reference frames file and either adds or remove a reference frame.
Current script is a non-specific addition/removal by selecting base pair frame at the mid-point of the configuration length.

Load a reference frames file, isolate a base pair frame in loop, and copy that right beside it.
Recall: For every n bp frame,
# - 1 = bp number and base-pair sequence [5*(n-1)+1]
# - 2 = origin position (x, y, z)        [5*(n-1)+2]
# - 3 = d1 vector                        [5*(n-1)+3]
# - 4 = d2 vector                        [5*(n-1)+4]
# - 5 = d3 vector                        [5*(n-1)+5]
'''
parser = argparse.ArgumentParser()
parser.add_argument("-in", "--input", dest="input", type=str, required=True, action='store', help="This arguement loads file that contains the list of temp configs.")
parser.add_argument("-fxn", "--function", choices=["increase", "decrease"], dest="function", type=str, required=True, action='store', help="Whether to increase or decrease the size of the DNA configuration")
parser.add_argument("-out", "--output", dest="outname", type=str, required=True, action='store', help="The name of the output files")
args = parser.parse_args()


infile = open( args.input, "r" )
lst = infile.readlines()
lst = [i.rstrip('\n') for i in lst]
infile.close()

if args.function=="increase":
    inheader     = lst[0]
    loop_len     = int(inheader.split()[0])
    midloop      = int( loop_len/2 )
    out_loop_len = int(loop_len+1)
    outheader    = " "+(str(out_loop_len))+" base pairs"
    new_loop = lst[ 1 : 5*(midloop)+1] + lst[5*(midloop-1)+1:5*(midloop)+1] + lst[5*(midloop)+1:]

elif args.function=="decrease":
    inheader     = lst[0]
    loop_len     = int(inheader.split()[0])
    midloop      = int( loop_len/2 )
    out_loop_len = int(loop_len - 1)
    outheader    = " "+(str(out_loop_len))+" base pairs"
    new_loop = lst[ 1 : 5*(midloop)+1] + lst[5*(midloop+1)+1:]


outfile = open(args.outname, "w")
outfile.write(outheader+'\n')
[outfile.write(i+'\n') for i in new_loop]
outfile.close()

del lst, inheader, loop_len, midloop, out_loop_len, outheader, new_loop
