import os
import re
import pyTestHarness.test as pthtest

# ---------------------------------------------------------------------------- #
def bottomUpSolid(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "bottomUpSolid"
  ranks = 1
  launch = [\
          'mkdir -p output',\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(thisDir,'test.opts'),\
          '../timeout.sh -t 4 ../../plot_spider.py 0',\
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
