#! /usr/bin/env python3

class extractInputParameter:
    "Pull a parameter out of the named file with parameterName(string), pPos, vPos as the name to look for, the position on the line where the name should occur and position on the line where the value should occur"


    def __init__(self):
        self.filename = "NoName"
        self.parameterName = "NoParameterName"
        self.vPos = -1
        self.pPos = -1


    def extract(self,filename,parameterName,vPos,pPos):
        self.filename = filename
        self.parameterName = parameterName
        self.vPos = vPos
        self.pPos = pPos

        parameterLines = list(open(self.filename))
        foundP = False
        value = 1e99
        for line in parameterLines:
            if self.parameterName in line:
                lineTokens = line.split()
                value = lineTokens[self.pPos]
                foundP = True
                continue

        if foundP == True:
            return value
        else:
            return "FAILED"


