import os
import re
import pyTestHarness.test as pthtest

# ---------------------------------------------------------------------------- #
def middleOutLiquid(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "middleOutLiquid"
  ranks = 1
  launch = [\
          'mkdir -p output',\
          os.path.join(rootDir,'main')  + ' -options_file ' + os.path.join(thisDir,'test.opts'),\
          '../../plot_figure.py 0,1000000',\
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
