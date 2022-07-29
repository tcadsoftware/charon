from __future__ import print_function

import math

class findCurrentAtParameter:
    "interpolate to find the current fiven a sweep parameter value"

#######################################################################################################
##  constructor
#######################################################################################################

    def __init__(self):
        self.parameters = []
        self.currents = []


#######################################################################################################
##  set data and check for sign
#######################################################################################################

    def setData(self,parameters,currents):
        for par in parameters:
            self.parameters.append(par)

        for lineNumber,c in enumerate(currents):
            self.currents.append(c)

#######################################################################################################
##  set target current
#######################################################################################################

    def setTargetParameter(self,targetParameter):
        self.targetParameter = targetParameter


#######################################################################################################
##  Interpolate to find the voltage
#######################################################################################################

    def calculateCurrent(self):
        #Find indexes that bracket the target
        # Note that this assumes paramters are monotonic
        startIndex = -1
        endIndex = -1

        descendingData = False
        if self.parameters[-1] < self.parameters[0]:
            descendingData = True

        for lineNumber,par in enumerate(self.parameters):
            if descendingData == False:
                if par >= self.targetParameter:
                    startIndex = lineNumber-1
                    endIndex = lineNumber
                    break
            if descendingData == True:
                if par <= self.targetParameter:
                    startIndex = lineNumber-1
                    endIndex = lineNumber
                    break



        #######################################################################################################
        # Now interpolate the voltage
        #######################################################################################################

        #If the denominator = zero, this will fail.
        if self.parameters[endIndex] == self.parameters[startIndex]:
            return (0,False)

        ratio = (self.targetParameter - self.parameters[startIndex])/(self.parameters[endIndex] - self.parameters[startIndex])

        current = self.currents[startIndex] + ratio*(self.currents[endIndex]-self.currents[startIndex])

        return (current,True)

    
