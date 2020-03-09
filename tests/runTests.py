#!/usr/bin/env python

import os
import sys
import traceback

thisDir = os.path.dirname(os.path.abspath(__file__))

# bitbucket.org/dmay/pythontestharness
sys.path.append(os.path.join(thisDir,'pythontestharness','lib'))  # overrides PYTHONPATH
try :
  import pyTestHarness.harness as pthharness
except Exception as errorMessage :
  if not sys.exc_info()[-1].tb_next :     # Check that the traceback has depth 1
    traceback.print_exc()
    print("********************")
    print("The required python library pyTestHarness was not found. Exiting.")
    print("If pyTestHarness is installed on your system, ensure pythontestharness/lib is included in the environment variable PYTHONPATH.")
    print("If pyTestHarness is not installed, obtain the source by executing the following:")
    print("  git clone https://bitbucket.org/dmay/pythontestharness " + os.path.join(thisDir,'pythontestharness'))
    print("********************")
    sys.exit(1)
  raise

# Main objects from pyTestHarness
import pyTestHarness.test as pthtest
import pyTestHarness.harness as pthharness
import pyTestHarness.version as pthversion

# Import functions to generate tests
from test_definitions import *

def main() :

    tol = 1.0E-9 # tolerance for checking values

    # Check that the version of pyTestHarness is great enough
    pthMajor,pthMinor,pthPath = pthversion.getVersion()
    if pthMinor < 3:
        raise RuntimeError("pyTestHarness version 0.3.0 or greater required")

    # Check for PETSc (should be the same as you used to compile)
    PETSC_DIR = os.getenv('PETSC_DIR')
    if not PETSC_DIR :
        print("You must define PETSC_DIR in your environment. Exiting.")
        sys.exit(1)
    PETSC_ARCH = os.getenv('PETSC_ARCH')
    if not PETSC_ARCH :
        print("You must define PETSC_ARCH in your environment. Exiting.")
        sys.exit(1)

    # The set of tests to run
    rootDir = os.path.join(thisDir,'..') # The SPIDER root directory

    allTests = [
        # the init tests below are not that useful, since they simply
        # return the initial condition. To make these useful requires
        # testing fields such as rho and temperature, rather than just
        # S and dS/dr
        #blackbody_init_liquid(rootDir, tol), \
        #blackbody_init_mixed(rootDir, tol), \
        #blackbody_init_solid(rootDir, tol), \
        # timestepping tests, below, are more useful
        blackbody(rootDir, tol) #,      \
        #atmosphere(rootDir, tol),     \
        #atmosphere_jeans(rootDir, tol)
    ]

    # Run tests
    os.environ['PYTHONUNBUFFERED'] = str('1')
    h = pthharness.Harness(allTests)
    h.execute()
    h.verify()

if __name__ == "__main__" :
    main()
