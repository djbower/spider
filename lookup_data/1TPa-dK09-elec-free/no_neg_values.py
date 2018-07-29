#!/usr/bin/env python

# script to ensure the data file has no negative entries in the first (pressure) column
# this zeros that entry and then writes the data to a new file.  Same number of decimal
# places is kept

# Eventually, the Jupityer notebooks of Aaron can do this operation, but for the time being,
# this gets us up and running

import sys
import numpy as np

infilename = sys.argv[1]

print( infilename )

infile = open( infilename, 'r' )
lines = infile.readlines()
infile.close()

outfile = open( infilename+'.zero', 'w' )

for nn, line in enumerate(lines):
    if nn<4:
        outfile.write( line )
    else:
        cols = line.split()
        col0 = cols[0]
        if float(col0) < 1.0E-6:
            col0n = '0.000000000000000000e+00'
            newline = '{} {} {}\n'.format( col0n, cols[1], cols[2] )
            outfile.write( newline )
        else:
            outfile.write( line )

outfile.close()

print( 'Done' )
