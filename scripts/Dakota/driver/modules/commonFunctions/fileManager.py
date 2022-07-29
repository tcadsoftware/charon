
import os
import sys

class fileManager:
    "file manager class opens, reads, closes"


    #######################################################################################################
    ##  file manager constructor
    #######################################################################################################

    def __init__(self):
        self.filename = "noName"
        self.inputLines = []
        self.vIndexFlag = False
        self.cIndexFlag = False


    #######################################################################################################
    ##  set filename to read
    #######################################################################################################

    def setFileName(self, filename):
        self.filename = filename
        # Check to make sure the file exists.
        if not os.path.isfile(self.filename):
            return False
        self.inputLines = list(open(self.filename))
        return True

    #######################################################################################################
    ##  set columns
    #######################################################################################################

    def setColumns(self, vIndex, voltString, cIndex, currentString):
        #For all Charon parameter sweeps, the first column will be varying voltage
        #self.voltColumn = 0
        self.vIndexFlag = False
        self.cIndexFlag = False
        if vIndex == "index":
            self.vIndexFlag = True
            self.voltColumn = int(voltString)

        if cIndex == "index":
            self.cIndexFlag = True
            self.currentColumn = int(currentString)



        #Split out the strings
        if self.vIndexFlag == False:
            voltList = voltString.split(",")

        if self.cIndexFlag == False:
            currentList = currentString.split(",")

        headerList = self.inputLines[0].split()
        # Note that the only time we have whitespace in a variable, it's in the parameter name
        # if the header has more columns than the data, it's because of the 0th column
        offset = len(headerList) - len(self.inputLines[1].split())
        if self.cIndexFlag == False:
            for index,header in enumerate(headerList):
                currentFound = True
                for cl in currentList:
                    if cl not in header.lower():
                        currentFound = False
                if currentFound == True:
                    self.currentColumn = index-offset

        if self.vIndexFlag == False:
            for index,header in enumerate(headerList):
                voltFound = True
                for vl in voltList:
                    if vl not in header.lower():
                        voltFound = False
                if voltFound == True:
                    self.voltColumn = index-offset


    #######################################################################################################
    ##  read file and store date in lists
    #######################################################################################################

    def readFile(self):
        volts = []
        current = []
        if self.voltColumn < 0 or self.currentColumn < 0:
            print ("Error:  Cannot determine current or voltage")
            sys.exit(1)

        #The first line of input is always text.
        for line in self.inputLines[1:]:
            lineTokens = line.split()
            volts.append(float(lineTokens[self.voltColumn]))
            current.append(float(lineTokens[self.currentColumn]))

        return (volts,current)


