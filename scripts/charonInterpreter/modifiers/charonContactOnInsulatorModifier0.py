
class charonContactOnInsulatorModifier0:
    "class for modifying the charonContactOnInsulatorModifier0 parameterList"


    def __init__(self):
        self.modifierName = "charonContactOnInsulatorModifier0"


    def getName(self):
        return self.modifierName


    def testForModification(self,pLList):
        foundMinValue = False
        foundMaxValue = False
        foundStepSize = False
        makeMinMaxModification = False
        makeStepSizeModification = False
        for lineNumber, line in enumerate(pLList):
            # Capture the min value
            if line.find("Charon->Solution Control->LOCA->Stepper,Min Value,double") >= 0:
                lineParts = line.split(",")
                minValue = lineParts[-1]
                minValueLine = lineNumber
                foundMinValue = True
            # Capture the max value
            if line.find("Charon->Solution Control->LOCA->Stepper,Max Value,double") >= 0:
                lineParts = line.split(",")
                maxValue = lineParts[-1]
                maxValueLine = lineNumber
                foundMaxValue = True
            # Capture the initial step size
            if line.find("Charon->Solution Control->LOCA->Step Size,Initial Step Size,double") >= 0:
                lineParts = line.split(",")
                initialStepSize = lineParts[-1]
                initialStepSizeLine = lineNumber
                foundStepSize = True
        if foundMinValue == True and foundMaxValue == True:
            if float(minValue) > float(maxValue):
                makeMinMaxModification = True
                replacementMinLine = "Charon->Solution Control->LOCA->Stepper,Min Value,double,"+maxValue
                replacementMaxLine = "Charon->Solution Control->LOCA->Stepper,Max Value,double,"+minValue
        if makeMinMaxModification and float(initialStepSize) > 0:
            makeStepSizeModification = True
            newStepSize = str(-float(initialStepSize))
            replacementStepSizeLine = "Charon->Solution Control->LOCA->Step Size,Initial Step Size,double,"+newStepSize
        # Make modifications
        if makeMinMaxModification == True:
            pLList[minValueLine] = replacementMinLine
            pLList[maxValueLine] = replacementMaxLine
        if makeStepSizeModification == True:
            pLList[initialStepSizeLine] = replacementStepSizeLine
        return pLList
