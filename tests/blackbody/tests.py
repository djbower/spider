import os
import re
import pyTestHarness.test as pthtest

def blackbody(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts'),\
          '../timeout.sh -t 6' + os.path.join(rootDir,'py3','plot_spider.py') + '-t 0,100,200,400,800,1200,1500',\
          ]

  # TODO stop depending on this
  expectedFile = os.path.join(thisDir,'output_expected','sol_100.m')

  # TODO
  def comparefunc(t) :
      t.compareFloatingPointRelative('',5e-1) # VERY loose tolerance!

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('output/sol_100.m') # TODO
  t.setVerifyMethod(comparefunc)
  t.appendKeywords('%')
  t.appendKeywords('[')
  t.appendKeywords(']')
  t.appendKeywords('sol')
  t.setWalltime(10) # minutes
  t.setUseSandbox()

  return(t)
