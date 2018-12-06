#! /usr/bin/env python
# load output JSON dump a txt file of a given timestep, for testing comparison purposes

from __future__ import print_function # maybe this will work with Python 2
import os
import sys
import json

timestep = 100

json_data = open(os.path.join('output',str(timestep)+'.json'))
data_d = json.load(json_data)
sol_values_array = data_d['solution']['values array']

outfile = open('out.txt','w')

for e in sol_values_array :
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

outfile.close()
json_data.close()
