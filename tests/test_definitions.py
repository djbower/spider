import os
import re
import pyTestHarness.test as pthtest

def atmosphere(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "atmosphere"
  ranks = 1
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(49990,50010)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','tests','atmosphere','bu_input.opts') + ' -nstepsmacro 1',\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_atmosphere.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '), tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(10) # minutes
  t.setUseSandbox()
  return(t)

def atmosphere_jeans(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "atmosphere_jeans"
  ranks = 1
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(49990,50010)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','tests','atmosphere_jeans','bu_input.opts') + ' -nstepsmacro 1',\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_atmosphere_jeans.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(10) # minutes
  t.setUseSandbox()
  return(t)

def blackbody(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','tests','blackbody','bu_input.opts'),\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py output/1049.json'), # TODO: last timestep time may change?
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_blackbody_1049.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(10) # minutes
  t.setUseSandbox()
  return(t)

# Currently not a useful test, since it only returns the initial condition
def blackbody_init_liquid(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody_init_liquid"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','tests','blackbody','bu_input.opts') + ' -ic_adiabat_entropy 2700' + ' -nstepsmacro 0',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_blackbody_init_liquid.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()

  return(t)

# Currently not a useful test, since it only returns the initial condition
def blackbody_init_mixed(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody_init_mixed"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','tests','blackbody','bu_input.opts') + ' -ic_adiabat_entropy 1900' + ' -nstepsmacro 0',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_blackbody_init_mixed.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()
  return(t)

# Currently not a useful test, since it only returns the initial condition
def blackbody_init_solid(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody_init_solid"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','tests','blackbody','bu_input.opts') + ' -ic_adiabat_entropy 1100' + ' -nstepsmacro 0',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_blackbody_init_solid.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()
  return(t)
