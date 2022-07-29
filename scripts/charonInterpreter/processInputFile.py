
import os
import sys
import subprocess
from os.path import exists,join

class processInputFile:
    "processes an input file for inlcudes, aprepro and maybe other stuff"

    def  __init__(self):
        self.inputFileLines = []
        self.processedIncludes = False


    def processIncludes(self,inputFileLines):
        self.inputFileLines = inputFileLines
        insertLocations = []
        insertTextFiles = []
        insertTextLines = []
        for index,line in enumerate(self.inputFileLines):
            lineTokens = line.split()
            if len(lineTokens) != 0:
                if lineTokens[0].lower() == "/include":
                    insertLocations.append(index)
                    insertTextFiles.append(lineTokens[1])

        offset = len(insertLocations)-1
        for index,insertion in enumerate(reversed(insertLocations)):
            #check for the existence of the inlcuded file
            rindex = offset-index
            if not os.path.isfile(insertTextFiles[rindex]):
                print ("Error!  I cannot find included file ",insertTextFiles[rindex])
                sys.exit(1)
            # pull in lines to add
            insertTextLines = list(open(insertTextFiles[rindex]))
            #remove the include line
            self.inputFileLines.pop(insertion)
            self.inputFileLines[insertion:insertion] = insertTextLines
            insertTextLines.clear()
        self.processedIncludes = True

        return self.inputFileLines

    def apreproInputFile(self,inputFileLines):
        #This seems silly, but we need to write out a temporary input file for aprepro
        # This assumes that the includes have already been processed.
        if not self.processedIncludes:
            print ("Error!! Cannot aprepro input file until includes have been processed")
            sys.exit(1)

        # First, write out a temp file for prepro
        inFilename = "temp.in"
        inFile = open(inFilename,'w+')
        for line in inputFileLines:
            inFile.write(line)
        inFile.close()

        # make sure we have aprepro
        apreproExecutable = "aprepro"
        haveAPREPRO = False
        searchPath = os.environ['PATH'].split(":")
        apreproPath = ""
        if 'CHARON_EXECUTABLE_PATH' in os.environ:
            searchPath.append(os.environ['CHARON_EXECUTABLE_PATH'])
        if 'APREPRO_PATH' in os.environ:
            searchPath.append(os.environ['APREPRO_PATH'])
        for sP in searchPath:
            if exists(join(sP,apreproExecutable)):
                haveAPREPRO = True
                apreproPath = sP+"/"

        if not haveAPREPRO:
            print ("aprepro not found.  Carrying on without preprocessing")
            if exists("temp.in"):
                os.remove("temp.in")
            return inputFileLines

        try:
            runCommand = apreproPath+"aprepro temp.in temp.out"
            returnValue = subprocess.check_call(runCommand,shell=True)

        except subprocess.CalledProcessError as e:
            print ("aprepro failed",e)

        returnInputLines = list(open("temp.out"))

        # Clean up files
        if exists("temp.in"):
            os.remove("temp.in")
        if exists("temp.out"):
            os.remove("temp.out")

        return returnInputLines

        
