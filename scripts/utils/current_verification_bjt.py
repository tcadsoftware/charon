###############################################################################
#
#  This script will attempt to verify the current and voltage values on the
#  base, collector, and emitter of a bipolar junction transistor.
#
#  The naming convention for the variables in this script is as follows:
#  * If it begins with 'i', it pertains to a current;
#    if it begins with 'v', it pertains to a voltage.
#  * If the second letter is 'B', it pertains to the base;
#    if the second letter is 'C', it pertains to the collector;
#    if the second letter is 'E', it pertains to the emitter.
#  * If it ends with 'Exact', it's an exact value specified by the user;
#    otherwise it's an estimated value computed by Charon.
#
###############################################################################

import getopt, sys, re
try:

  # Create a string to describe the usage of this script.
  usage = """\
Usage:  current_verification_bjt.py [inputs], where [inputs] are:
  -h|--help:              Display this usage information.
  -f|--filename:          The output file from a run of charon_mp.exe.
  -e|--emitterCurrent:    The true emitter current (in A/cm if in 2D, in A if
                          in 3D) against which to compare.
  -b|--baseCurrent:       The true base current (in A/cm if in 2D, in A if in
                          3D) against which to compare.
  -c|--collectorCurrent:  The true collector current (in A/cm if in 2D, in A if
                          in 3D) against which to compare.
  -t|--tolerance:         The tolerance for testing equality.
  --emitterVoltage:       The true voltage (in V) on the emitter if the current
                          constraint is placed there.  This input is optional.
  --baseVoltage:          The true voltage (in V) on the base if the current
                          constraint is placed there.  This input is optional.
  --collectorVoltage:     The true voltage (in V) on the collector if the
                          current constraint is placed there.  This input is
                          optional.
"""

  # Initialize the input variables.
  filename = ""
  iEExact = None
  iBExact = None
  iCExact = None
  tol = None
  vEExact = None
  vBExact = None
  vCExact = None

  # Get the command-line arguments.
  try:
    opts, args = getopt.getopt(
      sys.argv[1:], "hf:b:c:e:t:", [
        "help", "filename=", "emitterCurrent=", "baseCurrent=",
        "collectorCurrent=", "tolerance=", "emitterVoltage=", "baseVoltage=",
        "collectorVoltage="
      ]
    )
  except getopt.GetoptError:
    raise Exception(usage)
  for opt, arg in opts:
    if opt in ("-h", "--help"):
      print(usage)
      sys.exit()
    elif opt in ("-f", "--filename"):
      filename = arg
    elif opt in ("-e", "--emitterCurrent"):
      iEExact = float(arg)
    elif opt in ("-b", "--baseCurrent"):
      iBExact = float(arg)
    elif opt in ("-c", "--collectorCurrent"):
      iCExact = float(arg)
    elif opt in ("-t", "--tolerance"):
      tol = float(arg)
    elif opt == "--emitterVoltage":
      vEExact = float(arg)
    elif opt == "--baseVoltage":
      vBExact = float(arg)
    elif opt == "--collectorVoltage":
      vCExact = float(arg)

  # Ensure that we have all the command-line arguments we need.
  if filename == "":
    raise Exception("Must specify a filename.")
  if iEExact == None:
    raise Exception("Must specify an emitter current.")
  if iBExact == None:
    raise Exception("Must specify a base current.")
  if iCExact == None:
    raise Exception("Must specify a collector current.")
  if tol == None:
    raise Exception("Must specify a tolerance.")

  # Initialize the variables.
  iE  = None
  iB  = None
  iC  = None
  vE  = None
  vB  = None
  vC  = None
  dim = None

  # Create the current/voltage regular expressions for which to search.
  iEPat  = re.compile('    emitter_\S+_Current = (.*)$')
  iBPat  = re.compile('    base_\S+_Current = (.*)$')
  iCPat  = re.compile('    collector_\S+_Current = (.*)$')
  vEPat  = re.compile('  emitter\S+Voltage = (.*)$')
  vBPat  = re.compile('  base\S+Voltage = (.*)$')
  vCPat  = re.compile('  collector\S+Voltage = (.*)$')
  dimPat = re.compile('   Spatial dim = (.*)$')

  # Loop over the lines of the input file, trying to find the current/voltage
  # values.
  print("\nSearching the input file for currents/voltages...")
  f = open(filename)
  for line in f:
    iEMatch  = iEPat.match(line)
    iBMatch  = iBPat.match(line)
    iCMatch  = iCPat.match(line)
    vEMatch  = vEPat.match(line)
    vBMatch  = vBPat.match(line)
    vCMatch  = vCPat.match(line)
    dimMatch = dimPat.match(line)
    if iEMatch:
      iE = float(iEMatch.group(1))
    if iBMatch:
      iB = float(iBMatch.group(1))
    if iCMatch:
      iC = float(iCMatch.group(1))
    if vEMatch:
      vE = float(vEMatch.group(1))
    if vBMatch:
      vB = float(vBMatch.group(1))
    if vCMatch:
      vC = float(vCMatch.group(1))
    if dimMatch:
      dim = int(dimMatch.group(1))
  # end loop over the lines of the input file

  # Ensure that we found all the currents/voltages we need.
  if iE == None:
    raise Exception("Emitter current not found!")
  if iB == None:
    raise Exception("Base current not found!")
  if iC == None:
    raise Exception("Collector current not found!")
  if (vE == None) and (vEExact != None):
    raise Exception("Emitter voltage not found!")
  if (vB == None) and (vBExact != None):
    raise Exception("Base voltage not found!")
  if (vC == None) and (vCExact != None):
    raise Exception("Collector voltage not found!")

  # Print out the input arguments for debugging purposes.
  iUnit = ""
  if dim is 2:
    iUnit = " A/cm"
  elif dim is 3:
    iUnit = " A"
  print("Input Arguments:")
  print("File:                     " + filename)
  print("Exact Emitter Current:    " + str(iEExact) + iUnit)
  print("Exact Base Current:       " + str(iBExact) + iUnit)
  print("Exact Collector Current:  " + str(iCExact) + iUnit)
  print("Tolerance:                " + str(tol))
  if vEExact != None:
    print("Exact Emitter Voltage:    " + str(vEExact) + " V")
  if vBExact != None:
    print("Exact Base Voltage:       " + str(vBExact) + " V")
  if vCExact != None:
    print("Exact Collector Voltage:  " + str(vCExact) + " V")

  # Print out the results.
  print("\nResults:")
  print("Total Current     = " + str(iB + iC + iE).rjust(22) + iUnit)
  print("Emitter Current   = " + str(iE).rjust(22) + " (err = "
        + str(abs(iEExact - iE)).rjust(22) + ")" + iUnit)
  print("Base Current      = " + str(iB).rjust(22) + " (err = "
        + str(abs(iBExact - iB)).rjust(22) + ")" + iUnit)
  print("Collector Current = " + str(iC).rjust(22) + " (err = "
        + str(abs(iCExact - iC)).rjust(22) + ")" + iUnit)
  if vEExact != None:
    print("Emitter Voltage   = " + str(vE).rjust(22) + " (err = "
          + str(abs(vEExact - vE)).rjust(22) + ") V")
  if vBExact != None:
    print("Base Voltage      = " + str(vB).rjust(22) + " (err = "
          + str(abs(vBExact - vB)).rjust(22) + ") V")
  if vCExact != None:
    print("Collector Voltage = " + str(vC).rjust(22) + " (err = "
          + str(abs(vCExact - vC)).rjust(22) + ") V")

  # Ensure that the currents/voltages were within the specified tolerance.
  if abs(iEExact - iE) > tol:
    raise Exception("Emitter current value did not hit required tolerance.")
  if abs(iBExact - iB) > tol:
    raise Exception("Base current value did not hit required tolerance.")
  if abs(iCExact - iC) > tol:
    raise Exception("Collector current value did not hit required tolerance.")
  if vEExact != None: 
    if abs(vEExact - vE) > tol:
      raise Exception("Emitter voltage value did not hit required tolerance.")
  if vBExact != None: 
    if abs(vBExact - vB) > tol:
      raise Exception("Base voltage value did not hit required tolerance.")
  if vCExact != None: 
    if abs(vCExact - vC) > tol:
      raise Exception("Collector voltage value did not hit required "
                      + "tolerance.")

# If we throw any exceptions, that means the test has failed.  Otherwise it's
# passed.
except Exception as e:
  print("Test Failed: ")
  print(e)
else:
  print("Test Passed")
