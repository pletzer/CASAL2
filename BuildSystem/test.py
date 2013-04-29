import os
import sys
import os.path
import subprocess
import fileinput

sys.path.insert(0, "buildtools/classes")

# These are iSAM specific python objects
from System import *
from Globals import *
from Builder import *

"""
Print the usage for this build system
"""
def print_usage():  
  os.system( [ 'clear', 'cls' ][ os.name == 'nt' ] )
  print '###########################################################'
  print '# iSAM Build System Usage                                 #'
  print '###########################################################'
  print 'Usage:'
  print 'doBuild <build_target> <build_parameter>'
  print ''
  print 'Valid Build Types:'
  print '  debug - Build debug version'
  print '  release - Build release version'
  print '  thirdparty - Build all required third party libraries'
  print '  all - Does clean, thirdparty, debug, release builds in order'
  print '  clean - Remove any previous debug/release build information'
  print '  cleanall - Remove all previous build information'
  print ''
  print 'Valid Build Parameters: (thirdparty only)'
  print '  <libary name> - Target third party library to build or rebuild'
  print ''
  print 'Valid Build parameters: (debug/release only)'
  print '  admb - Use ADMB auto-differentiation in compiled executable'


"""
Start our build system. This method essentially just
gets all of the system information we need and will
test that executables are in the path that we need
"""
def start_build_system():
  system_info.set_new_path()
  if not system_info.find_compiler_path():
    return False
  
  #if Globals.operating_system_ == "win32":
    #system_info.find_cmd_path()
    
  return True  

"""
Start The main execution here:
"""
system_info = SystemInfo()
if not start_build_system():
  system_info.reset_original_path()
  exit(1)

"""
Get the build information from the user
"""
def start():
  build_target = ""
  build_parameters = ""
  
  allowed_build_targets = [ "debug", "release", "thirdparty", "test", "all", "clean", "cleanall", "help" ]
 
  if not len(sys.argv) > 1:  
    os.system( [ 'clear', 'cls' ][ os.name == 'nt' ] )
    print "*************************************************************************"
    print "*************************************************************************"
    print "Welcome to the iSAM Build System"
    print "*************************************************************************"
    print "*************************************************************************"
    print "\n"
    print "Valid build targets: " +  ", ".join(allowed_build_targets)
    build_target = raw_input("Please enter a valid build type [debug]: ").lower()
    while build_target != "" and not build_target in allowed_build_targets:
      print "## Error: " + build_target + " is not a valid build target" 
      build_target = raw_input("Please enter a valid build type [debug]: ").lower()
    if build_target == "":
      build_target = "debug" 
    if build_target == "help":
      print_usage()
      return True
  
    # Find our build parameters to use, this will enable admb etc
    if build_target == "debug" or build_target == "release":    
      print "Valid build parameters: " + ", ".join(allowed_build_parameters)
      build_parameters = raw_input("Please enter build parameters [none]: ").lower()
      while build_parameters != "" and not build_parameters in allowed_build_parameters:
        print "## Error: " + build_parameters + " is not a valid build parameter"
        build_parameters = raw_input("Please enter build parameters [none]: ").lower()
      if build_parameters == "none":
        build_parameters = ""  
  
  """
  Handle build information already passed in
  """
  if len(sys.argv) > 1 and len(str(sys.argv[1])) > 1:
      build_target = sys.argv[1]
  if len(sys.argv) > 2 and len(str(sys.argv[2])) > 1:
      build_parameters = sys.argv[2] 
   
  if build_target == "" or not build_target.lower() in allowed_build_targets:
    return Globals.PrintError(build_target + " is not a valid build target")
    
  build_target = build_target.lower()    
  if build_target == "help":
    print_usage()
    return True

  if build_parameters != "": 
    build_parameters = build_parameters.lower()
  
  Globals.build_target_ = build_target
  Globals.build_parameters_ = build_parameters
  
  print " -- Build target: " + Globals.build_target_
  print " -- Build parameters: " + Globals.build_parameters_
  
  if build_target == "debug" or build_target == "release":
    allowed_build_parameters = [ "", "admb" ]    
    if not build_parameters in allowed_build_parameters:
      return Globals.PrintError("Build parameter " + build_parameters + " is not valid")
    
    os.system( [ 'clear', 'cls' ][ os.name == 'nt' ] )
    print "*************************************************************************"
    print "*************************************************************************"
    print "--> Starting " + Globals.build_target_ + " Build"
    code_builder = MainCode()
    if not code_builder.start():
      return False
  elif build_target == "thirdparty":
    os.system( [ 'clear', 'cls' ][ os.name == 'nt' ] )
    print "*************************************************************************"
    print "*************************************************************************"
    print "--> Starting " + Globals.build_target_ + " Build"
    code_builder = ThirdPartyLibraries()
    if not code_builder.start():
      return False      
      
"""
This is the entry point in to our build system
"""
exit_code = 0
if not start():
  exit_code = 1
  
system_info.reset_original_path()
print "--> Finished"
exit(exit_code)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 

  