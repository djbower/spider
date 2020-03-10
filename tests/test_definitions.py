import os
import re
import pyTestHarness.test as pthtest

def atmosphere(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "atmosphere"
  ranks = 1
  #acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(49990,50010)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','atmosphere.opts') + ' -nstepsmacro 1',\
          #os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py output/50001.json'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_atmosphere_50001.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '), tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('50001.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(2) # minutes
  t.setUseSandbox()
  return(t)

def atmosphere_ic(rootDir, tol) :
   thisDir = os.path.split(os.path.abspath(__file__))[0]
   testName = "atmosphere_ic"
   ranks = 1
   launch = [\
           os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','atmosphere.opts') + ' -nstepsmacro 0',\
           os.path.join(rootDir,'tests','json_timestep_to_txt.py output/0.json'),
            ]

   expectedFile = os.path.join(thisDir,'expected_output','expected_atmosphere_0.txt')

   def comparefunc(t) :
       t.compareFloatingPointRelative(re.escape('scaling: '), tol)
       t.compareFloatingPointRelative(re.escape('val: '), tol)

   t = pthtest.Test(testName,ranks,launch,expectedFile)
   t.setComparisonFile('0.txt')
   t.setVerifyMethod(comparefunc)
   t.setWalltime(1) # minutes
   t.setUseSandbox()
   return(t)

# TODO: escape test
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
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','blackbody.opts'),\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py output/1005.json'), # TODO: last timestep time may change?
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_blackbody_1005.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('1005.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()
  return(t)

def reaction(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "reaction"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','reaction.opts') + ' -nstepsmacro 1',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py output/10000.json'), # TODO: last timestep time may change?
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_reaction_10000.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('10000.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(2) # minutes
  t.setUseSandbox()
  return(t)

def reaction_ic(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "reaction_ic"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','reaction.opts') + ' -nstepsmacro 0',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py output/0.json'), # TODO: last timestep time may change?
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_reaction_0.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('0.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()
  return(t)
