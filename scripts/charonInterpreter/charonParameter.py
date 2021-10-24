from __future__ import print_function

class charonParameter:
    "parameter can be string, double, int, bool"

    def __init__(self):
        self.parameterName = "Value"
        self.parameterType = "string"
        self.parameterValue = ""
        self.satisfied = "False"

    def printParameter(self,indention):
        self.printLine = indention+"<Parameter name=\""+self.parameterName+"\" type=\""+self.parameterType+"\" value=\""+self.parameterValue+"\" />"
        print(self.printLine,end='\n')

    def getParameterString(self,indention):
        self.printLine = indention+"<Parameter name=\""+self.parameterName+"\" type=\""+self.parameterType+"\" value=\""+self.parameterValue+"\" />"
        return self.printLine

    def addMe(self,parameterName,parameterType,parameterValue):
        self.parameterName  = parameterName
        self.parameterType  = parameterType
        self.parameterValue = parameterValue
        self.satisfied = "True"
        return self

    def replaceParameterValue(self,parameterType,parameterValue):
        self.parameterType  = parameterType
        self.parameterValue = parameterValue
        self.satisfied = "True"

    def errorCheckParameter(self):
        return (self.satisfied,self.parameterName)


