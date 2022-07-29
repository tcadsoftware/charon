
import sys

from modules.commonFunctions.parse.parse import *
from parseExecuteMethod import *
from processInputFile import *

class preProcessCommand:
    "process a list from the execute methods to do preprocessing prior to execution"


    def __init__(self):
        self.executeMethod = []
        self.preProcessCommand = []
        self.pI = processInputFile()

    def preProcess(self,inputParametersFilename,executeMethod):
        # This is simple
        self.preProcessCommand.clear()
        self.executeMethod = executeMethod
        self.inputParametersFilename = inputParametersFilename
        app = ""
        template = ""
        options = ""
        pEM = parseExecuteMethod()

        (app,template,options) = pEM.parseEM(self.executeMethod)



        # in every case, a template should be provided.  The suffix of the 
        # template will be template, but the active file produced may have various 
        # suffixes.  Can check for those.

        filenameSuffix = ""

        if app.lower() == "charon":
            filenameSuffix = "inp"

        if app.lower() == "cubit" or app.lower() == "pymesh":
            filenameSuffix = "jou"

        for temp in template:
            #Process any includes and create a temporary working file for dprepro
            tempInput = list(open(temp))
            processedTemplateLines =self.pI.processIncludes(tempInput)
            processedTemplatesString = " ".join(processedTemplateLines)
            processedTemplatesFilename = temp+"consolidated"
            pTF = open(processedTemplatesFilename,'w+')
            pTF.write(processedTemplatesString)
            pTF.close()
            self.preProcessCommand.append("dprepro "+inputParametersFilename+" "+processedTemplatesFilename+" "+temp.replace("template",filenameSuffix))

        return self.preProcessCommand


