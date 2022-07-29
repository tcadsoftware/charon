
import sys
import os
from os import path
from os import listdir
from os.path import exists, join
import subprocess

class domainDecomposition:
    "Provide for domain decomposition and re-decomposition"

    def __init__(self):
        self.filename = ""
        self.foundDecomp = False
        self.foundEpu = False
        self.pathToDecomp = ""
        self.newSize = 1
        self.oldSize = 1
        #Get the tools
        (self.decompTool,self.epuTool) = self.searchForTools()
        self.decompTool += "/decomp"
        self.epuTool += "/epu"


    def searchForTools(self):
        searchPath = os.environ['PATH'].split(":")
        decomp = ""
        epu = ""
        for sP in searchPath:
            if exists(join(sP,"decomp")):
                decomp = sP
                self.foundDecomp = True
            if exists(join(sP,"epu")):
                epu = sP
                self.foundEpu = True
        return (decomp,epu)


    def decompExists(self):
        return self.foundDecomp

    def epuExists(self):
        return self.foundEpu


    def decompose(self,filename,numProcs):
        self.filename = filename
        self.newSize = numProcs
        # Check if decomposition exists
        decompFile = filename+"."+str(numProcs)+".0"
        if exists(decompFile):
            print("Decomposition of this size already exists.\n Please remove it before proceeding.")
            return False

        try:
            runCommand = self.decompTool
            arguments = " --processors "+str(self.newSize)+" "+self.filename
            totalCommand = runCommand+arguments
            print ("executing:  ",runCommand+arguments)
            returnValue = subprocess.check_call(totalCommand,shell=True)
            return True

        except subprocess.CalledProcessError as e:
            print ("Decomposition failed.  Aborting.",e)
            return False

    def resize(self,filename,numProcs,oldNumProcs):
        self.filename = filename
        self.newSize = numProcs
        self.oldSize = oldNumProcs
        # Check if decomposition of old size exists--This is a little bit goofy way to do this, but other things would be more cumbersome
        decompFile = filename+"."+str(oldNumProcs)+".0"
        decompFile1 = decompFile+"0"
        decompFile2 = decompFile1+"0"
        decompFile3 = decompFile2+"0"
        if not exists(decompFile) and not exists(decompFile1) and not exists(decompFile2) and not exists(decompFile3):
            print("Cannot find ",decompFile,"Please make certain the options are correct.")
            return False

        # Check if decomposition of new size exists--This is a little bit goofy way to do this, but other things would be more cumbersome
        decompFile = filename+"."+str(numProcs)+".0"
        decompFile1 = decompFile+"0"
        decompFile2 = decompFile1+"0"
        decompFile3 = decompFile2+"0"
        if exists(decompFile) or exists(decompFile1) or exists(decompFile2) or exists(decompFile3):
            print("Decomposition of size ",str(numProcs)+" already exists.\n Please remove it before proceeding.")
            return False

        try:
            runCommand = self.epuTool
            arguments = " -auto "+self.filename+"."+str(self.oldSize)+"."
            totalCommand = runCommand+arguments
            print ("executing:  ",runCommand+arguments)
            returnValue = subprocess.check_call(totalCommand,shell=True)

        except subprocess.CalledProcessError as e:
            print ("epu'ing failed.  Aborting.",e)
            return False

        try:
            runCommand = self.decompTool
            arguments = " --processors "+str(self.newSize)+" "+self.filename
            totalCommand = runCommand+arguments
            print ("executing:  ",runCommand+arguments)
            returnValue = subprocess.check_call(totalCommand,shell=True)
            return True

        except subprocess.CalledProcessError as e:
            print ("Decomposition failed.  Aborting.",e)
            return False



