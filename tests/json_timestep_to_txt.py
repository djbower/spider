#! /usr/bin/env python
# load output JSON dump a txt file of a given timestep, for testing comparison purposes

from __future__ import print_function # maybe this will work with Python 2
import os
import sys
import json

# Interpret the first argument as a timestep to look for
if len(sys.argv) < 2 :
    timestep = 0
else :
    timestep = int(sys.argv[1])

with open(os.path.join('output',str(timestep)+'.json')) as json_data :
    data_d = json.load(json_data)
    subdomain_data_array = data_d['solution']['subdomain data']
    with open('out.txt','w') as outfile :
        for e in subdomain_data_array :
            outfile.write('description: ')
            outfile.write(e['description'])
            outfile.write('\n')
            outfile.write('scaling: ')
            outfile.write(e['scaling'])
            outfile.write('\n')
            for ee in e['values'] :
                outfile.write('val: ')
                outfile.write(ee)
                outfile.write('\n')
