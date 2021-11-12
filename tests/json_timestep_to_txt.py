#! /usr/bin/env python
# load output JSON dump a txt file of a given timestep, for testing comparison purposes

from __future__ import print_function # maybe this will work with Python 2
import os
import sys
import json

# Interpret the first argument as a JSON file to look for
if len(sys.argv) < 2 :
    filename = 'output/0.json'
else :
    filename = sys.argv[1]

timestep = filename.split('/')[-1].split('.json')[0]
outfilename = '{}.txt'.format(timestep)

with open(filename) as json_data:
    data_d = json.load(json_data)
    subdomain_data_array = data_d['solution']['subdomain data']
    with open(outfilename,'w') as outfile :
        for e in subdomain_data_array:
            # only compare values that exist
            # for example, if atmosphere and reactions are turned off then
            # some values are not available to compare 
            val_l = e['values']
            if len(val_l) == 0:
                continue
            outfile.write('description: ')
            outfile.write(e['description'])
            outfile.write('\n')
            outfile.write('scaling: ')
            outfile.write(e['scaling'])
            outfile.write('\n')

            for nn, ee in enumerate(e['values']) :
                outfile.write('val: ')
                outfile.write(ee)
                outfile.write('\n')
