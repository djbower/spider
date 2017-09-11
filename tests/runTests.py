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

# Import functions to generate tests
from bottomUpLiquid.tests  import *
from bottomUpMixed.tests   import *
from bottomUpSolid.tests   import *
from middleOutLiquid.tests import *
from middleOutMixed.tests  import *
from middleOutSolid.tests  import *

# ---------------------------------------------------------------------------- #
def main() :

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
    rootDir = os.path.join(thisDir,'..')
    allTests = [
        bottomUpLiquid(rootDir), \
        bottomUpMixed(rootDir),  \
        bottomUpSolid(rootDir),  \
        middleOutLiquid(rootDir),\
        middleOutMixed(rootDir), \
        middleOutSolid(rootDir), \
    ]

    # Run tests
    os.environ['PYTHONUNBUFFERED'] = str('1')
    h = pthharness.Harness(allTests)
    h.execute()
    h.verify()

# ---------------------------------------------------------------------------- #
if __name__ == "__main__" :
    main()
