#! /usr/bin/env python3

import math
import re

class pyMeshPrepro:
    "Do prepro arithmetic on refinement directives"


    def __init__(self,cubitObj):
        self.cubitObj = cubitObj


    def processLine(self,line):
        variableList = list(self.cubitObj.get_aprepro_vars())
        returnString = line

        if "{" in line:
            apVar = re.findall(r'{(.+?)}',line)
            apVarStr = " ".join(apVar)
            apVarStrOriginal = apVarStr
            # Process the line by inserting spaces around arithmetic operators
            apVarStr = " ".join(apVarStr.replace("+"," + ").split())
            apVarStr = " ".join(apVarStr.replace("-"," - ").split())
            apVarStr = " ".join(apVarStr.replace("*"," * ").split())
            apVarStr = " ".join(apVarStr.replace("/"," / ").split())

            # replace the apVar List
            apVar = apVarStr.split()

            for index,ap in enumerate(apVar):
                if ap in variableList:
                    apVar[index] = self.cubitObj.get_aprepro_numeric_value(ap)
                elif ap == "+":
                    apVar[index] = self.addition
                elif ap == "-":
                    apVar [index]= self.subtraction
                elif ap == "*":
                    apVar[index] = self.multiplication
                elif ap == "/":
                    apVar[index] = self.division
                else:  #assume its a double if nothing else
                    apVar[index] = float(ap)

            returnValue = apVar[0]
            for idx,ap in enumerate(apVar[1:]):
                index = idx+1
                if type(ap) != float:
                    returnValue = ap(returnValue,apVar[index+1])

            #Now make the processed replacements in the directive line
            returnString = returnString.replace("{","")
            returnString = returnString.replace("}","")
            returnString = returnString.replace(apVarStrOriginal,str(returnValue))

        return returnString

    def addition(self,val1,val2):
        return val1+val2


    def subtraction(self,val1,val2):
        return val1-val2

    def multiplication(self,val1,val2):
        return val1*val2

    def division(self,val1,val2):
        return val1/val2
