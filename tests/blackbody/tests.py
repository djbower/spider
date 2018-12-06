import os
import re
import pyTestHarness.test as pthtest

def blackbody(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody"
  ranks = 1
  # Run the test,
  # then use a local script to generate a text file to compare with,
  # then run the plotting script (output not checked)
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts'),\
          os.path.join(thisDir,'json_timestep_to_txt.py'),
          os.path.join(rootDir,'tests','timeout.sh') + '-t 6' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0,100,200,400,800,1200,1500',\
          ]

  expectedFile = os.path.join(thisDir,'expected.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '),1e-14)
      t.compareFloatingPointRelative(re.escape('val: '),    1e-10) # likely too tight

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(10) # minutes
  t.setUseSandbox()

  return(t)
