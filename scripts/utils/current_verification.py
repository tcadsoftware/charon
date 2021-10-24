###############################################################################
#
#  This script will attempt to verify the current and voltage values on the
#  anode and cathode of a diode.
#
#  The naming convention for the variables in this script is as follows:
#  * If it begins with 'i', it pertains to a current;
#    if it begins with 'v', it pertains to a voltage.
#  * If the second letter is 'A', it pertains to the anode;
#    if the second letter is 'C', it pertains to the cathode.
#  * If it ends with 'Exact', it's an exact value specified by the user;
#    otherwise it's an estimated value computed by Charon.
#
###############################################################################

import getopt, sys, re
try:

  # Create a string to describe the usage of this script.
  usage = """\
Usage:  current_verification.py [inputs], where [inputs] are:
  -h|--help:            Display this usage information.
  -f|--filename:        The output file from a run of charon_mp.exe.
  -a|--anodeCurrent:    The true anode current (in A/cm if in 2D, in A if in
                        3D) against which to compare.
  -c|--cathodeCurrent:  The true cathode current (in A/cm if in 2D, in A if in
                        3D) against which to compare.
  -t|--tolerance:       The tolerance for testing equality.
  --anodeVoltage:       The true voltage (in V) on the anode if the current
                        constraint is placed there.  This input is optional.
  --cathodeVoltage:     The true voltage (in V) on the cathode if the current
                        constraint is placed there.  This input is optional.
"""

  # Initialize the input variables.
  filename = ""
  iAExact = None
  iCExact = None
  tol = None
  vAExact = None
  vCExact = None

  # Get the command-line arguments.
  try:
    opts, args = getopt.getopt(
      sys.argv[1:], "hf:a:c:t:", [
        "help", "filename=", "anodeCurrent=", "cathodeCurrent=", "tolerance=",
        "anodeVoltage=", "cathodeVoltage="
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
    elif opt in ("-a", "--anodeCurrent"):
      iAExact = float(arg)
    elif opt in ("-c", "--cathodeCurrent"):
      iCExact = float(arg)
    elif opt in ("-t", "--tolerance"):
      tol = float(arg)
    elif opt == "--anodeVoltage":
      vAExact = float(arg)
    elif opt == "--cathodeVoltage":
      vCExact = float(arg)

  # Ensure that we have all the command-line arguments we need.
  if filename == "":
    raise Exception("Must specify a filename.")
  if iAExact == None:
    raise Exception("Must specify an anode current.")
  if iCExact == None:
    raise Exception("Must specify a cathode current.")
  if tol == None:
    raise Exception("Must specify a tolerance.")

  # Initialize the variables.
  iA  = None
  iC  = None
  vA  = None
  vC  = None
  dim = None

  # Create the regular expressions for which to search.
  iAPat  = re.compile('    anode_\S+_Current = (.*)$')
  iCPat  = re.compile('    cathode_\S+_Current = (.*)$')
  vAPat  = re.compile('  anode\S+Voltage = (.*)$')
  vCPat  = re.compile('  cathode\S+Voltage = (.*)$')
  dimPat = re.compile('   Spatial dim = (.*)$')

  # Loop over the lines of the input file, trying to find the current/voltage
  # values.
  print("\nSearching the input file for currents/voltages...")
  f = open(filename)
  for line in f:
    iAMatch  = iAPat.match(line)
    iCMatch  = iCPat.match(line)
    vAMatch  = vAPat.match(line)
    vCMatch  = vCPat.match(line)
    dimMatch = dimPat.match(line)
    if iAMatch:
      iA = float(iAMatch.group(1))
    if iCMatch:
      iC = float(iCMatch.group(1))
    if vAMatch:
      vA = float(vAMatch.group(1))
    if vCMatch:
      vC = float(vCMatch.group(1))
    if dimMatch:
      dim = int(dimMatch.group(1))
  # end loop over the lines of the input file

  # Ensure that we found all the currents/voltages we need.
  if iA == None:
    raise Exception("Anode current not found!")
  if iC == None:
    raise Exception("Cathode current not found!")
  if (vA == None) and (vAExact != None):
    raise Exception("Anode voltage not found!")
  if (vC == None) and (vCExact != None):
    raise Exception("Cathode voltage not found!")

  # Print out the input arguments for debugging purposes.
  iUnit = ""
  if dim is 2:
    iUnit = " A/cm"
  elif dim is 3:
    iUnit = " A"
  print("Input Arguments:")
  print("File:                   " + filename)
  print("Exact Anode Current:    " + str(iAExact) + iUnit)
  print("Exact Cathode Current:  " + str(iCExact) + iUnit)
  print("Tolerance:              " + str(tol))
  if vAExact != None:
    print("Exact Anode Voltage:    " + str(vAExact) + " V")
  if vCExact != None:
    print("Exact Cathode Voltage:  " + str(vCExact) + " V")

  # Print out the results.
  print("\nResults:")
  print("Total Current   = " + str(iA + iC).rjust(22) + iUnit)
  print("Anode Current   = " + str(iA).rjust(22) + " (err = "
        + str(abs(iAExact - iA)).rjust(22) + ")" + iUnit)
  print("Cathode Current = " + str(iC).rjust(22) + " (err = "
        + str(abs(iCExact - iC)).rjust(22) + ")" + iUnit)
  if vAExact != None:
    print("Anode Voltage   = " + str(vA).rjust(22) + " (err = "
          + str(abs(vAExact - vA)).rjust(22) + ") V")
  if vCExact != None:
    print("Cathode Voltage = " + str(vC).rjust(22) + " (err = "
          + str(abs(vCExact - vC)).rjust(22) + ") V")

  # Ensure that the currents/voltages were within the specified tolerance.
  if abs(iAExact - iA) > tol:
    raise Exception("Anode current value did not hit required tolerance.")
  if abs(iCExact - iC) > tol:
    raise Exception("Cathode current value did not hit required tolerance.")
  if vAExact != None: 
    if abs(vAExact - vA) > tol:
      raise Exception("Anode voltage value did not hit required tolerance.")
  if vCExact != None: 
    if abs(vCExact - vC) > tol:
      raise Exception("Cathode voltage value did not hit required tolerance.")

# If we throw any exceptions, that means the test has failed.  Otherwise it's
# passed.
except Exception as e:
  print("Test Failed: ")
  print(e)
else:
  print("Test Passed")
