
import math

class findParameterAtLogCurrent:
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

    def setData(self,parameters,current):
        for v in parameters:
            self.parameters.append(v)

        for lineNumber,c in enumerate(current):
            c = math.fabs(c)
            clog = math.log10(c)
            self.current.append(clog)

#######################################################################################################
##  set target current
#######################################################################################################

    def setTargetCurrent(self,targetCurrent):
        self.targetCurrent = targetCurrent
        self.targetCurrentLog = math.log10(targetCurrent)



#######################################################################################################
##  Interpolate to find the parameter
#######################################################################################################

    def calculateParameter(self):
        #Find indexes that bracket the target
        startIndex = -1
        endIndex = -1

        for lineNumber,c in enumerate(self.current):
            if c > self.targetCurrentLog:
                startIndex = lineNumber-1
                endIndex = lineNumber
                break

        if self.current[endIndex] == self.current[startIndex]:
            return (0.0,False)

        #######################################################################################################
        # Now interpolate the parameter
        #######################################################################################################
        ratio = (self.targetCurrentLog - self.current[startIndex])/(self.current[endIndex] - self.current[startIndex])

        parameter = self.parameters[startIndex] + ratio*(self.parameters[endIndex]-self.parameters[startIndex])

        return (parameter,True)

    
