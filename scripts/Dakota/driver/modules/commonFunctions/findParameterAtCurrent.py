
import math

class findParameterAtCurrent:
    "interpolate to find the parameter at a given current"

#######################################################################################################
##  constructor
#######################################################################################################

    def __init__(self):
        self.parameters = []
        self.current = []


#######################################################################################################
##  set data and check for sign
#######################################################################################################

    def setData(self,parameters,current,useAbs=False):
        self.parameters.extend(parameters)

        if useAbs == True:
            for lineNumber,c in enumerate(current):
                self.current.append(math.fabs(c))
        else:
            self.current.extend(current)

#######################################################################################################
##  set target current
#######################################################################################################

    def setTargetCurrent(self,targetCurrent):
        self.targetCurrent = targetCurrent

#######################################################################################################
##  Interpolate to find the parameter
#######################################################################################################

    def calculateParameter(self):
        #Find indexes that bracket the target
        startIndex = -1
        endIndex = -1

        for lineNumber,c in enumerate(self.current):
            if c > self.targetCurrent:
                startIndex = lineNumber-1
                endIndex = lineNumber
                break

        if self.current[endIndex] == self.current[startIndex]:
            return (0.0,False)

        #######################################################################################################
        # Now interpolate the parameter
        #######################################################################################################
        ratio = (self.targetCurren - self.current[startIndex])/(self.current[endIndex] - self.current[startIndex])

        parameter = self.parameters[startIndex] + ratio*(self.parameters[endIndex]-self.parameters[startIndex])

        return (parameter,True)

    
