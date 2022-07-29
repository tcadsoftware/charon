# When the user disables the WRMS test in a run using tpetra the
# associated status test needs to be modified to reflect the fact that
# all that is left is a NormF test.
#
# Note you can't remove lines in the list at this point because the
# interpreter makes assumptions about the number of lines being fixed
# by the time modifiers are invoked. Instead I'm setting the lines to
# the empty string.
start Modifier 0

def testForModification(self, pList):

  # If the nonlinear solver tolerance was specified by the user pull
  # it out here and save it
  foundRythmos = False
  for entry in pList:
    if entry.find("Charon->Solution Control->NOX->Status Tests->Test 0->Test 0,Tolerance,double,") != -1:
      foundRythmos = True
      tolLineRythmos = entry

  # There is possibly a similar line for tempus so handle that as well
  foundTempus  = False
  for entry in pList:
    if entry.find("Charon->Solution Control->Tempus->Default Stepper->Default Solver->NOX->Status Tests->Test 0->Test 0,Tolerance,double,") != -1:
      foundTempus = True
      tolLineTempus = entry

  # Tolerance should actually always be found since it's set to a
  # default value if it isn't specified, but test to avoid future
  # issues.
  #
  # Pull out the value of Tolerance and put it into the base "Test 0"
  # instead of "Test 0->Test 0".
  tolVal = float(1.0e-4)
  if foundRythmos:
    toks = tolLineRythmos.split(',')
    tolVal = float(toks[len(toks)-1])

  # Get rid of the sub tests in Test 0 since Test 0 will be a scalar
  # NormF test when WRMS is disabled.
  pList[:] = [entry if entry.find("Charon->Solution Control->NOX->Status Tests->Test 0->") == -1 else str("") for entry in pList]

  # Replace the entries left over from when WRMS was on
  pList[:] = [entry if entry.find("Charon->Solution Control->NOX->Status Tests->Test 0,Test Type,string,Combo") == -1 else str("") for entry in pList]
  pList[:] = [entry if entry.find("Charon->Solution Control->Tempus->Default Stepper->Default Solver->NOX->Status Tests->Test 0,Test Type,string,Combo") == -1 else str("") for entry in pList]
  
  pList[:] = [entry if entry.find("Charon->Solution Control->NOX->Status Tests->Test 0,Combo Type,string,AND") == -1 else str("") for entry in pList]
  pList[:] = [entry if entry.find("Charon->Solution Control->Tempus->Default Stepper->Default Solver->NOX->Status Tests->Test 0,Combo Type,string,AND") == -1 else str("") for entry in pList]

  pList[:] = [entry if entry.find("Charon->Solution Control->NOX->Status Tests->Test 0,Number of Tests,int,2") == -1 else str("") for entry in pList]
  pList[:] = [entry if entry.find("Charon->Solution Control->Tempus->Default Stepper->Default Solver->NOX->Status Tests->Test 0,Number of Tests,int,2") == -1 else str("") for entry in pList]

  # Now we need to replace a line just blanked with the valid
  # Tolerance specification
  foundRythmos = False
  for i in range(len(pList)):
    if pList[i] == "":
      pList[i] = str("Charon->Solution Control->NOX->Status Tests->Test 0,Tolerance,double,")+str(tolVal)
      pList[i] = str("Charon->Solution Control->Tempus->Default Stepper->Default Solver->NOX->Status Tests->Test 0,Tolerance,double,")+str(tolVal)
      foundRythmos = True
      break

  return pList

end Modifier 0
