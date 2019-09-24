import os
import re
import pyTestHarness.test as pthtest

def atmosphere(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "atmosphere"
  ranks = 1
  # Long range of acceptable files (probably better to check earlier)
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(1000000,1000124)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','atmosphere','bu_input.opts'),\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'),
          #os.path.join(rootDir,'tests','timeout.sh') + ' -t 6 ' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0,100,200,400,800,1200,1500',\
          ]

  expectedFile = os.path.join(thisDir,'expected_atmosphere.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '),1e-14)
      t.compareFloatingPointRelative(re.escape('val: '),    1e-5)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(10) # minutes
  t.setUseSandbox()
  return(t)

def atmosphere_jeans(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "atmosphere_jeans"
  ranks = 1
  # Long range of acceptable files (probably better to check earlier)
  acceptable_files = ' '.join(['output/'+str(i)+'.json' for i in range(1000000,1000097)])
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','atmosphere_jeans','bu_input.opts'),\
          os.path.join(rootDir,'tests','copy_one_of.py') + ' ' + acceptable_files + ' -o to_check.json',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py to_check.json'),
          #os.path.join(rootDir,'tests','timeout.sh') + ' -t 6 ' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0,100,200,400,800,1200,1500',\
          ]

  expectedFile = os.path.join(thisDir,'expected_atmosphere_jeans.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '),1e-14)
      t.compareFloatingPointRelative(re.escape('val: '),    1e-5)

  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setComparisonFile('out.txt')
  t.setVerifyMethod(comparefunc)
  t.setWalltime(10) # minutes
  t.setUseSandbox()
  return(t)

def blackbody(rootDir) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = "blackbody"
  ranks = 1
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts'),\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py output/100.json'),
          #os.path.join(rootDir,'tests','timeout.sh') + ' -t 6 ' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0,100,200,400,800,1200,1500',\
          ]

  expectedFile = os.path.join(thisDir,'expected.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '),1e-14)
      t.compareFloatingPointRelative(re.escape('val: '),    1e-4)

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
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts') + ' -ic_adiabat_entropy 2700' + ' -nstepsmacro 0',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py'),
          #os.path.join(rootDir,'tests','timeout.sh') + ' -t 6 ' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0',\
          ]

  expectedFile = os.path.join(thisDir,'expected_init_liquid.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '),1e-14)
      t.compareFloatingPointRelative(re.escape('val: '),    1e-5)

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
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts') + ' -ic_adiabat_entropy 1900' + ' -nstepsmacro 0',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py'),
          #os.path.join(rootDir,'tests','timeout.sh') + ' -t 6 ' + os.path.join(rootDir,'py3','plot_spider.py') + ' -t 0',\
          ]

  expectedFile = os.path.join(thisDir,'expected_init_mixed.txt')

  def comparefunc(t) :
      t.compareFloatingPointRelative(re.escape('scaling: '),1e-14)
      t.compareFloatingPointRelative(re.escape('val: '),    1e-5)

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
  launch = [\
          os.path.join(rootDir,'spider')  + ' -options_file ' + os.path.join(rootDir,'examples','bower_2019','blackbody','bu_input.opts') + ' -ic_adiabat_entropy 1100' + ' -nstepsmacro 0',\
          os.path.join(rootDir,'tests','json_timestep_to_txt.py'),
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
