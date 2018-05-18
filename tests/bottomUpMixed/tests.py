import os
import re
import pyTestHarness.test as pthtest

# ---------------------------------------------------------------------------- #
def bottomUpMixed(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "bottomUpMixed"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(thisDir,'test.opts'),\
          '../timeout.sh -t 6 ../../py3/plot_bower_et_al_2018.py -t 0 -f3',\
          ]

  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
      # TODO
      pass

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setVerifyMethod(comparefunc)
  t.setUseSandbox()

  return(t)
