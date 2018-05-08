import os
import re
import pyTestHarness.test as pthtest

def bottomUpLiquid(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "bottomUpLiquid"
  ranks = 1
  launch = [\
          'mkdir -p output',\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(thisDir,'test.opts'),\
          '../timeout.sh -t 6 ../../py3/plot_spider.py 0,60,120,180,240 -f3',\
          ]

  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
      t.compareFloatingPointRelative('',5e-1) # VERY loose tolerance!

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('output/sol_20.m')
  t.setVerifyMethod(comparefunc)
  t.appendKeywords('%')
  t.appendKeywords('[')
  t.appendKeywords(']')
  t.appendKeywords('sol')
  t.setUseSandbox()

  return(t)
