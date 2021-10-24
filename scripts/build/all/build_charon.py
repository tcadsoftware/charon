#!/usr/bin/env python3

"""
Requires Python >= 2.4

Building Charon on my (GLH) workstation (Ubuntu 12.04, gnu 4.4.x
compilers), and on my Macbook air (Darwin 11.4.1 (lion), gnu 4.4.x
compilers)

Current as of 18Dec2012

Prerequisites:

 (read Panzer/drekar/build_scripts/README_BuildingTPLS.txt for
  further explanation)

 1) MPI (required)
 2) Boost (required)
 3) netcdf (required)
 4) hdf5 (optional)
 5) zlib (optional)
 6) ParMetis (optional)
 7) PAPI (optional)
 8) Drekar/Panzer cmake flags to correctly enable the above support

To Build:

 1) check out Trilinos using:
        git clone software.sandia.gov:/space/git/Trilinos

 2) cd into the Trilinos directory and check out charon:
        cd Trilinos
        git clone software.sandia.gov:/space/sandiagit/charon

 3) create a build directory

 4) cd into the build directory and invoke this script

 5) invoke make

 6) Profit! (or not)
"""

#from __future__ import print_function

import sys
import os
import re
from optparse import OptionParser, OptionGroup
import platform
import subprocess
import getpass

###########################################################################
class BuildOpts:

  option_list = []
  keyvals = []
  resetvals = []

  ########################################################
  def __init__(self, dbglvl, keyvals, resetvals):
    self.dbglvl = dbglvl
    self.keyvals = keyvals
    self.resetvals = resetvals

  ########################################################
  def readFromFile(self, optString):
    option_list.add(optString)
    return

  ########################################################
  def writeToFile(self, fname):
    f = open(fname, 'w')
    for opt in self.option_list:
      if opt[1] == None:
        opt_str = opt[0] + '=' + opt[2] + '\n'
      else:
        opt_str = opt[0] + ':' + opt[1] + '=' + opt[2] + '\n'

      f.write(opt_str)

    f.close

    return

  ########################################################
  def addOption(self, optstring):

    tup = []

    # Only split based on the first equal symbol. This avoids
    # problems when the option itself has an equal symbol in it.
    line_vals = re.split('=', optstring, 1)

    if self.dbglvl >= 9:
      print("[DBG]: " + str(line_vals), len(line_vals))

    if len(line_vals) == 2: # var=val, no type

      # Now determine if it has a type specification
      line_sub_vals = re.split(':', line_vals[0])

      if len(line_sub_vals) == 1:
        tup = [line_vals[0], None, line_vals[1]]
      elif len(line_sub_vals) == 2:
        tup = [line_sub_vals[0], line_sub_vals[1], line_vals[1]]
      else:
        print("ERROR. Unhandled option found. Argument was\n\t" + optstring,
              file=sys.stderr)
        sys.exit(3)

    elif opstring.strip() != "":
      # anything else is an error
      err_msg = 'Unknown option format found\n'
      err_msg += fname + '\n'
      err_msg += 'Option was "' + str(optstring) + '"'
      raise Exception(err_msg)

    tup = self._check_and_replace_kv(tup)
    tup = self._check_and_replace_rv(tup)
    self.option_list.append(tup)

    return

  ########################################################
  def readFromFile(self, fname):
    file = open(fname)

    tup = []
    for line in file:
      line=line.strip('\n')

      # If the non-empty line starts with a '#' then it's a comment
      # and will be ignored
      if len(line) > 0 and line[0] != '#':
        self.addOption(line)

    file.close()

    if self.dbglvl >= 8:
      print('[DBG]: Options read were:')
      for tup in self.option_list:
        print('  [DBG]' + str(tup))

    return

  ########################################################
  def _check_and_replace_kv(self, in_array):

    out_array = in_array
    for keyval in self.keyvals:
      fcount = 0
      for field in in_array:
        if field:
          find_str = '\$\{'+keyval[0]+'\}'
          out_array[fcount] = re.sub(find_str, keyval[1], field)
        fcount = fcount + 1

    return out_array

  ########################################################
  def _check_and_replace_rv(self, in_array):

    if len(self.resetvals) > 0:
      for keyval in self.resetvals:
        if re.search(keyval[0], in_array[0]):
          in_array[2] = keyval[1]
          # Set a flag in the structure to indicate it was used so
          # that an error check can be performed later
          keyval[2] = '1'

    return in_array

  ########################################################
  def generateCMakeInvocation(self, cmake_cmd, cmake_opts, src_dir):

    # Find a cmake binary
    bin_path = os.environ['PATH']
    for dirname in bin_path.split(os.pathsep):
      fname = os.path.join(dirname, cmake_cmd)
      if os.path.exists(fname):
        found = True
        break
      else:
        found = False

    if not found:
      raise Exception('Could not find cmake in your path')

    if cmake_opts == None:
      cmake_opts = ''

    cmake_cmd_line = fname
    for opt in self.option_list:
      if opt[1] == None:
        add_str = ' -D ' + opt[0] + '=' + opt[2]
      else:
        add_str = ' -D ' + opt[0] + ':' + opt[1] + '=' + opt[2]

      cmake_cmd_line += add_str

    cmake_cmd_line += ' ' + str(cmake_opts) + ' ' + str(src_dir)

    return cmake_cmd_line;

###########################################################################
def main():

  usage = "usage: %prog [options] <tcad-charon directory>"
  parser = OptionParser(usage=usage)

  parser.add_option("-m", "--manual-page", action="store_true", dest="mpage",
                    help="print extended help")

  parser.add_option("-f", "--options-file", action="append",
                    dest="opfilenames",
                    help="file, or files, containing cmake options for build",
                    metavar="<existing file name>")

  parser.add_option("-b", "--build-type", type="choice",
                    action="store", dest="bldtype",
                    choices=["r", "release", "d", "debug"],
                    default="r",
                    help="type of build, debug or release")

  parser.add_option("--cmake-options", dest="cmakeopts",
                    help="options to pass directly to cmake invocation",
                    metavar="'STRING'")

  parser.add_option("-r", "--reset-option", action="append", dest="resetopts",
                    help="option to reset", metavar="SUBSTRING:NEW VALUE")

  parser.add_option("-d", "--dump-options-file", dest="dumpfilename",
                    help="file to store resulting cmake options for later use",
                    metavar="<file name to be created or overwritten>")

  parser.add_option("-k", "--key-value-replace", action="append",
                    dest="keyvals",
                    help="key, value replacement for templated options",
                    metavar="KEY:VALUE")

  parser.add_option("-c", "--options-directory", dest="optdir",
                    help="directory with option files",
                    metavar="<directory containing *.opts files>")

  parser.add_option("--cmake-command", dest="cmakecmd",
                    help="full path spec to cmake executable",
                    default="cmake",
                    metavar="<cmake command>")

  group = OptionGroup(parser, "Debugging options",
                      "These are for use when debugging this script")

  group.add_option("--echo-cmake-only", dest="echocmakecmd", default=False,
                    action="store_true",
                    help="just echo the cmake command; don't run cmake")

  group.add_option("--debug-level", dest="dbglvl", type=int, default=0,
                    help="level of debugging for this script",
                    metavar="<0-9>")

  parser.add_option_group(group)

  (options, args) = parser.parse_args()

  # Default relative path where this script resides
  script_rel_path = "/scripts/build/all"

  # Debug output of command line options
  if options.dbglvl >= 5:
    print("[DBG]: Options dump...")

    if options.opfilenames:
      for fname in options.opfilenames:
        print("    [DBG]: " + str(fname))
    else:
      print("  [DBG]: No option files given")

    if options.dumpfilename:
      print("  [DBG]: dump file name: " + str(options.dumpfilename))

    if options.resetopts:
      print("  [DBG]: options to reset: ")
      for optrst in options.resetopts:
        print("    [DBG]: " + str(optrst))

    if options.keyvals:
      print("  [DBG]: substring, replacement:")
      for kval in options.keyvals:
        print("    [DBG]: " + str(kval))

  # If the user asked for the manual page then print it and exit
  if options.mpage:
    print_man_page()
    sys.exit(0)

  # Check and parse, when necessary, the specified options
  opfilenames = []

  # Get the user name to check for a user-specific script
  # directory. This will be the new default if the user hasn't
  # overridden it via the command line.
  loginid = getpass.getuser()

  if options.dbglvl >= 3:
    print("[DBG]: loginid = " + str(loginid))
    print("[DBG]: sys.path[0] = " + sys.path[0])

  # Look for the base directory of personal opts files. This assumes
  # that the build_charon.py script is being run from the directory
  # containing the distribution that is going to be built,
  # specifically scripts/build/all.
  m = re.search(r'(.*?)\/all', sys.path[0])
  if m and loginid:
    optsfiledir = m.group(1) + '/' + loginid
  else:
    optsfiledir = '.'

  # Override the default options file directory if the user gave a
  # specific directory
  if options.optdir:
    optsfiledir = options.optdir

  if options.dbglvl >= 3:
    print("[DBG]: options file directory: " + optsfiledir)

  # If there is a "General.opts" file in the current directory it will
  # override everything else and any other *.opts files will be
  # ignored. This allows a user to use a single file with all their
  # options in it, similar to what was done with configure scripts.
  if os.path.exists("General.opts"):
    opfilenames.append("General.opts")
  else:
    # Make sure the path is valid
    if not os.path.exists(optsfiledir):
      print("ERROR: Specified directory\n   " + optsfiledir +
            "\nfor *.opts files does not exist", file=sys.stderr)
      sys.exit(2)

    # Add the opts files. Note that at least one file must be found
    if os.path.exists(optsfiledir+"/General.opts"):
      opfilenames.append(optsfiledir+"/General.opts")

    # Check for a platform specific file directory and a machine
    # specific file in the user directory
    test_fn = optsfiledir + '/' + platform.system() + '.opts'
    if options.dbglvl >= 1:
      print("[DBG]: Platform-specific file name: " + test_fn)
    if os.path.exists(test_fn):
      opfilenames.append(test_fn)

    tmp_host = platform.uname()[1]
    # Get rid of any domain name. Just want unqualified name
    hostname = tmp_host.split('.')[0]
    test_fn = optsfiledir + '/' + hostname + '.opts'
    if options.dbglvl >= 1:
      print("[DBG]: host-specific file name: " + test_fn)
    if os.path.exists(test_fn):
      opfilenames.append(test_fn)

    if len(opfilenames) == 0:
      print("ERROR: Could not find any *.opts files in the " +
            "specified directory\n   " + optsfiledir, file=sys.stderr)
      sys.exit(2)

    # Check for the environment variable SNLSYSTEM and use that for an
    # opts file if it exists.
    snl_sys = os.environ.get('SNLSYSTEM')

    if snl_sys:
      if options.dbglvl >= 1:
        print("[DBG]: value of SNLSYSTEM: " + snl_sys)

    if snl_sys:
      test_fn = optsfiledir + '/' + snl_sys + '.opts'
      if options.dbglvl >= 1:
        print("[DBG]: SNL specific file name: " + test_fn)
      if os.path.exists(test_fn):
        opfilenames.append(test_fn)

    # Now add any user specified option files to the list. These will
    # take precedence over the others. Any options specified within
    # will override anything also found in the preceeding files
    if options.opfilenames:
      for fname in options.opfilenames:
        if not os.path.exists(fname):

          # Try prepending the user directory
          if not os.path.exists(optsfiledir+"/"+fname):
            print("ERROR: Could not find options file names " + fname,
              file=sys.stderr)
            sys.exit(3)
          else:
            opfilenames.append(optsfiledir+"/"+fname)
        else:
          opfilenames.append(fname)


  if options.dbglvl >= 2:
    if opfilenames:
      print("[DBG]: Found option files...")
      for opfile in opfilenames:
        print("  [DBG]: " + opfile)

  # If the script doesn't find any option files report that as an error
  if not opfilenames:
    print("ERROR: No opt files found. You must specify at least one " +
          "using the -f option", file=sys.stderr)
    sys.exit(2)


  if options.dbglvl >= 1:
    print ('[DBG]: args = ', args)

  if not args:
    tcad_charon_dir = re.sub(script_rel_path, '', sys.path[0])
  elif len(args) > 1:
    print("ERROR:  The only argument allowed is the the path to " +
          "the tcad-charon directory to utilize for the build.", file=sys.stderr)
    sys.exit(4)
  else:
    tcad_charon_dir = args[0]

  if options.dbglvl >= 1:
    print ('[DBG]: tcad_charon_dir = ', tcad_charon_dir)

  # Put the key-value pairs into a nested list
  kv_list = []
  if options.keyvals:
    for kv in options.keyvals:
      kv_list.append(re.split(':', kv))

  if options.dbglvl >= 8:
    print('[DBG]: kv_list = ', kv_list)

  # Put reset options into a nested list
  rv_list = []
  if options.resetopts:
    for rv in options.resetopts:
      rv = rv + ":0"
      rv_list.append(re.split(':', rv))

  if options.dbglvl >= 8:
    print('[DBG]: rv_list = ', rv_list)

  # Only positional argument allowed is the tcad-charon directory
  if not os.path.exists(tcad_charon_dir):
    print("ERROR:  The tcad-charon directory \"" + str(tcad_charon_dir) +
          "\" doesn't exist.", file=sys.stderr)
    sys.exit(5)

  # Everything should be good to go at this point. Note that the order
  # of this loop determines the priority. In the most general case the
  # priority should be 1) <hostname>.opts, 2) <system>.opts and 3)
  # General.opts in the user directory. General.opts in the current
  # directory means none of the others are used.
  bld_opts = BuildOpts(options.dbglvl, kv_list, rv_list)
  for fname in opfilenames:
    if options.dbglvl >= 5:
      print('[DBG]: reading options file: ' + fname)

    bld_opts.readFromFile(fname)

  # Add the build type cmake option
  if options.bldtype == 'r' or options.bldtype == 'release':
    bld_opts.addOption("CMAKE_BUILD_TYPE:STRING=RELEASE")
  else:
    bld_opts.addOption("CMAKE_BUILD_TYPE:STRING=DEBUG")

  if options.dumpfilename:
    bld_opts.writeToFile(options.dumpfilename)

  # Check to make sure all requested resets were performed
  if len(rv_list) > 0:
    for keyval in rv_list:
      if keyval[2] != '1':
        raise Exception("ERROR: The requested option, " + keyval[0]
                        + ", for reset was not found")

  cmake_cmd = bld_opts.generateCMakeInvocation(options.cmakecmd,
                                               options.cmakeopts,
                                               tcad_charon_dir)

  # Remove the cache file and directory
  rm_cmd = '/bin/rm -f CMakeCache.txt'
  ret_code = subprocess.call(rm_cmd, shell=True)

  if options.dbglvl >= 9:
    print('[DBG]: cmake command would be...\n' + cmake_cmd)
    ret_code = 0
  elif options.echocmakecmd == True:
    print(cmake_cmd)
    ret_code = 0
  else:
    ret_code = subprocess.call(cmake_cmd, shell=True)

  sys.exit(ret_code)

###########################################################################
def print_man_page():

  print(" ")
  print("This script is used to invoke cmake for building Charon.")
  print(" ")
  print("The idea of this script is to have a set of files with the required cmake build")
  print("options. One file, General.opts, just contains general options you'll need for")
  print("any machine. The <system>.opts file contains system specific options where")
  print("<system> is an operating system, Linux or Darwin, for example. And lastly")
  print("a <host>.opts file containing options specific to the host named <host>.")
  print(" ")
  print("File searching:")
  print("  The script looks for files in the following order")
  print("    (1) A file named General.opts in the current directory")
  print("          NOTE: If this file exists no other files will be searched for.")
  print("    (2) A file named General.opts in the user directory")
  print("    (3) A file named <system>.opts in the user directory")
  print("    (4) A file named <hostname>.opts in the user directory")
  print("    (5) A file named <snlsystem>.opts in the user directory")
  print("    (6) Any script file, or files, specified by the user.")
  print("        where:")
  print("           <system> is the name of your OS, like Darwin (Mac OSX)")
  print("            or Linux, for example.")
  print("           <hostname> is the unqualified network name of the host.")
  print("           <snlsytem> comes from the environment variable SNLSYSTEM if it exists.")
  print(" ")
  print("     The \"user directory\" is currently given by \"tcad-charon/scripts/build/<login>\"")
  print("     <login> is your login ID")
  print(" ")
  print("     By default the files are assumed to be in the directory in which the script resides,")
  print("     however this can be overriden with the options --options-directory argument.")
  print(" ")
  print("     If a variable is set in multiple files then the last variable read will override ")
  print("     any previous variables. In this wey you can reset variables specified in previous ")
  print("     by redefining them in a later *.opts file")
  print(" ")
  print("Key-value replacement (-k or --key-value option):")
  print("  The specified options files may have keys delimited by ${<KEY>} where <KEY> is some")
  print("  key to be replaced by the value as specified on the command line. For example if the")
  print("  option \"-k NETCDFDIR:/usr/local\" were specified all the options would be searched")
  print("  for ${NETCDFDIR} and that would be replaced with \"/usr/local\"")
  print(" ")
  print("Resetting existing options (-r or --reset-option):")
  print("  To reset an existing option to a new value. For example, if the option")
  print("  tcad-charon_VERBOSE_CONFIGURE:BOOL=ON is specified you can reset it off with")
  print("  -r VERBOSE_CONFIGURE:OFF. Note that a substring search is used so any option with")
  print("  \"VERBOSE_CONFIGURE\" would be reset to off")
  print(" ")

  return

###########################################################################
if __name__ == "__main__":

  # Check version. Need at least 2.4 for "subprocess"
  ver = sys.version_info
  if ver < (2,4,0):
    print("ERROR: This script requires python >= 2.4, found version " +
          str(ver[0]) + "." + str(ver[1]), file=sys.stderr)
    sys.exit(1)

  main()
