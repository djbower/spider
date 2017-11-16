import os
import re
import pyTestHarness.test as pthtest

# ---------------------------------------------------------------------------- #
def bottomUpMixed(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "bottomUpMixed"
  ranks = 1
  launch = [\
          'mkdir -p output',\
          os.path.join(rootDir,'main')  + ' -options_file ' + os.path.join(thisDir,'test.opts'),\
          '../../plot_spider.py 0',\
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
