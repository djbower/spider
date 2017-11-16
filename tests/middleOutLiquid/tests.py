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
          '../../plot_spider.py 0,1000000,2000000',\
          ]

  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
      t.compareFloatingPoint('',1e-5) # match everything

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('output/dSdr_b_aug_1000000.m')
  t.setVerifyMethod(comparefunc)
  t.appendKeywords('%')
  t.appendKeywords('[')
  t.appendKeywords(']')
  t.appendKeywords('dSdr_b_aug')
  t.setUseSandbox()

  return(t)
