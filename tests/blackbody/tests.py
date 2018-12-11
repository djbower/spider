import os
import re
import pyTestHarness.test as pthtest

# Note: there is a lot of code duplication here. This should be addressed
# before copying this same logic any more

def blackbody(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody"
  ranks = 1
  # Run the test,
  # then use a local script to generate a text file to compare with,
  # then run the plotting script (output not checked)
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts'),\
          os.path.join(thisDir,'json_timestep_to_txt.py 100'),
          #os.path.join(rootDir,'tests','timeout.sh') + ' -t 6 ' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0,100,200,400,800,1200,1500',\
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

def blackbody_init_liquid(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody_init_liquid"
  ranks = 1
  # Run the test,
  # then use a local script to generate a text file to compare with,
  # then run the plotting script (output not checked)
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts') + ' -ic_adiabat_entropy 2700' + ' -nstepsmacro 0',\
          os.path.join(thisDir,'json_timestep_to_txt.py'),
          #os.path.join(rootDir,'tests','timeout.sh') + ' -t 6 ' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0',\
          ]

  expectedFile = os.path.join(thisDir,'expected_init_liquid.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '),1e-14)
      t.compareFloatingPointRelative(re.escape('val: '),    1e-10) # likely too tight

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()

  return(t)

def blackbody_init_mixed(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody_init_mixed"
  ranks = 1
  # Run the test,
  # then use a local script to generate a text file to compare with,
  # then run the plotting script (output not checked)
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts') + ' -ic_adiabat_entropy 1900' + ' -nstepsmacro 0',\
          os.path.join(thisDir,'json_timestep_to_txt.py'),
          #os.path.join(rootDir,'tests','timeout.sh') + ' -t 6 ' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0',\
          ]

  expectedFile = os.path.join(thisDir,'expected_init_mixed.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '),1e-14)
      t.compareFloatingPointRelative(re.escape('val: '),    1e-10) # likely too tight

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()

  return(t)

def blackbody_init_solid(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody_init_solid"
  ranks = 1
  # Run the test,
  # then use a local script to generate a text file to compare with,
  # then run the plotting script (output not checked)
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts') + ' -ic_adiabat_entropy 1100' + ' -nstepsmacro 0',\
          os.path.join(thisDir,'json_timestep_to_txt.py'),
          #os.path.join(rootDir,'tests','timeout.sh') + ' -t 6 ' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0',\
          ]

  expectedFile = os.path.join(thisDir,'expected_init_solid.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '),1e-14)
      t.compareFloatingPointRelative(re.escape('val: '),    1e-10) # likely too tight

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()

  return(t)
