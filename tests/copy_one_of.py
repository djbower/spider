#!/usr/bin/env python
from __future__ import print_function
import os
import shutil
import argparse

def copy_one_of(filenames_in,filename_out) :
    filename_found = None
    for filename_in in filenames_in :
        if os.path.isfile(filename_in) :
            if filename_found :
                raise Exception("More than one input file exists.")
            filename_found = filename_in
    if not filename_found :
        raise Exception("No input files found")
    shutil.copyfile(filename_found,filename_out)

def main() :
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename-out','-o',required=True)
    args,filenames_in = parser.parse_known_args()
    copy_one_of(filenames_in,args.filename_out)

if __name__ == "__main__" :
    main()
