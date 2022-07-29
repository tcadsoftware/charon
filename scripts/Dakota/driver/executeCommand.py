
import sys
import os

from modules.commonFunctions.parse.parse import *
from parseExecuteMethod import *

class executeCommand:
    "Create an execute command based on execute methods"


    def __init__(self):
        self.parametersDict = {}
        self.executeMethod = []
        if 'CHARON_EXECUTABLE_PATH' in os.environ:
            self.executablePath = os.environ['CHARON_EXECUTABLE_PATH']+"/"
        else:
            self.executablePath = ""

    def execute(self,parametersDict,executeMethod):
        self.parametersDict = parametersDict
        self.executeMethod = executeMethod
        pEM = parseExecuteMethod()

        (app,template,options) = pEM.parseEM(self.executeMethod)

        executor = "None"
        inputFilename = "noName"

        #Handle the execute method (cubit, pyMesh, chirp)
        if app.lower() == "cubit":
            executor = "cubit"
            inputFilename = template[0].replace("template","jou")
        elif app.lower() == "pymesh":
            executor = self.executablePath + "pyMesh"
            inputFilename = template[0].replace("template","jou")
        elif app.lower() == "charon":
            executor = self.executablePath + "chirp -i"
            inputFilename = template[0].replace("template","inp")
        else:
            print ("Error: No known executor found for ",app)
            sys.exit(1)
            
        executeCommand = executor+" "+inputFilename

        runLabel = ""
        if 'label' in parametersDict:
            runLabel = " -l evalID-"+parametersDict['label']+" "

        #Add extras for chirp
        if app.lower() == "charon":
            executeCommand += " --silent --run --cleanTextData "+runLabel+" --np "+str(self.parametersDict['executeProcs'])

        #Add extras for chirp
        if app.lower() == "pymesh":
            executeCommand += " |& tee "+inputFilename+".log >& /dev/null"

        #Append the execute command with any arguments that may have been supplied
        executeCommand += " "+options

        #Return the final command
        return executeCommand

