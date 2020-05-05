import os
import re
import pyTestHarness.test as pthtest

def atmosphere(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "atmosphere"
  ranks = 1
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(49950,50050)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','atmosphere.opts') + ' -nstepsmacro 1',\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_atmosphere.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '), tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('to_check.txt')
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

   expectedFile = os.path.join(thisDir,'expected_output','expected_atmosphere_ic.txt')

   def comparefunc(t) :
       t.compareFloatingPointRelative(re.escape('scaling: '), tol)
       t.compareFloatingPointRelative(re.escape('val: '), tol)

   t = pthtest.Test(testName,ranks,launch,expectedFile)
   t.setComparisonFile('0.txt')
   t.setVerifyMethod(comparefunc)
   t.setWalltime(1) # minutes
   t.setUseSandbox()
   return(t)

def atmosphere_escape(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "atmosphere_escape"
  ranks = 1
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(49950,50050)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','atmosphere_escape.opts') + ' -nstepsmacro 1',\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_atmosphere_escape.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('to_check.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(2) # minutes
  t.setUseSandbox()
  return(t)

def atmosphere_escape_ic(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "atmosphere_escape_ic"
  ranks = 1
  #acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(49990,50010)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','atmosphere_escape.opts') + ' -nstepsmacro 0',\
          #os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py output/0.json'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_atmosphere_escape_ic.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('0.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()
  return(t)

def atmosphere_escape_zerosol(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "atmosphere_escape_zerosol"
  ranks = 1
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(49750,50250)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','atmosphere_escape_zerosol.opts') + ' -nstepsmacro 1',\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'),
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_atmosphere_escape_zerosol.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('to_check.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(2) # minutes
  t.setUseSandbox()
  return(t)

def blackbody(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody"
  ranks = 1
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(950,1050)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','blackbody.opts'),\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'), # TODO: last timestep time may change?
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_blackbody.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('to_check.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()
  return(t)

def reaction(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "reaction"
  ranks = 1
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(9050,10050)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','reaction.opts') + ' -nstepsmacro 1',\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'), # TODO: last timestep time may change?
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_reaction.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('to_check.txt')
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

  expectedFile = os.path.join(thisDir,'expected_output','expected_reaction_ic.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('0.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()
  return(t)

def reaction_zerosol(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "reaction_zerosol"
  ranks = 1
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(9050,10050)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','reaction_zerosol.opts') + ' -nstepsmacro 1',\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'), # TODO: last timestep time may change?
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_reaction_zerosol.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('to_check.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(2) # minutes
  t.setUseSandbox()
  return(t)

def reaction_zerosol_ic(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "reaction_zerosol_ic"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','reaction_zerosol.opts') + ' -nstepsmacro 0',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py output/0.json'), # TODO: last timestep time may change?
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_reaction_zerosol_ic.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('0.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()
  return(t)

def socrates_p_restart(rootDir, tol) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "socrates_p_restart"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'tests','opts','socrates_p_restart.opts')
              + ' -ic_interior_filename ' + os.path.join(rootDir,'tests','opts','socrates_p_restart_63745.json'),\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py output/90611.json'), # TODO: last timestep time may change?
          ]

  expectedFile = os.path.join(thisDir,'expected_output','expected_socrates_p_restart.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '), tol)
      t.compareFloatingPointRelative(re.escape('val: '),    tol)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('90611.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(1) # minutes
  t.setUseSandbox()
  return(t)
