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
            # DJB: no need to compare values if they do not exist
            # for example, if atmosphere and reactions are turned off then
            # some values are not solved for 
            val_l = e['values']
            if len(val_l) == 0:
                continue
            # DJB: description could change, but shouldn't change the numerical result
            outfile.write('description: ')
            outfile.write(e['description'])
            outfile.write('\n')
            # DJB: scaling could change, but shouldn't change the (physical) numerical result
            outfile.write('scaling: ')
            outfile.write(e['scaling'])
            outfile.write('\n')
            for ee in e['values'] :
                outfile.write('val: ')
                outfile.write(ee)
                outfile.write('\n')
