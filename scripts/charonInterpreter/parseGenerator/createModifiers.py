
import os


class createModifiers:
    "create modifiers"

#######################################################################################################
##  Create modifiers from the full line of modifiers in the parseInputs directory structure
#######################################################################################################
 
    def __init__(self,modifierSourcePath,modifierTargetPath,verbosity):
        self.modifierSourcePath = modifierSourcePath
        self.modifierTargetPath = modifierTargetPath
        self.verbosity = verbosity
        self.modifierLibraryList = []



#######################################################################################################
##  Loop over directory structure and
#######################################################################################################

    def generateModifiers(self):
        modifierInputList = []
        modifierInputFilenameList = []

        # Collect up all filenames
        for dirname, dirnames, filenames in os.walk(self.modifierSourcePath):
            for filename in filenames:
                if filename[-14:].lower() == "modifierinp.py":
                    modifierInputList.append(os.path.join(dirname,filename))
                    modifierInputFilenameList.append(filename)

        # Iterate over the modifier list and create modifiers
        for modNumber, mIFL in enumerate(modifierInputFilenameList):
            modifierCanonicalName = "charon"+mIFL.replace(mIFL[-15:],"").replace(mIFL[:6],"")
            self.modifierLines = list(open(modifierInputList[modNumber]))
            modifierFunctions = self.extractFunctions()
            #write modifiers
            for modNumber, mF in enumerate(modifierFunctions):
                modifierName = modifierCanonicalName+"Modifier"+str(modNumber)
                self.modifierLibraryList.append(modifierName)
                modifierFunction = self.getCompletedClass(modifierName,mF)
                modifierFile = open(self.modifierTargetPath+"/"+modifierName+".py","w+")
                modifierFile.write(modifierFunction)
                modifierFile.close()


#######################################################################################################
##  extract modifier functions
#######################################################################################################

    def extractFunctions(self):

        functionCounter = -1
        functionOpen = False
        functionLists = []

        for line in self.modifierLines:
            lineTokens = line.split()

            if len(lineTokens) == 0:
                continue

            if line[0] == "#":
                continue

            if functionOpen == True:
                if lineTokens[0].lower() != "end":
                    functionLists[functionCounter].append(line)
                if lineTokens[0].lower() == "end":
                    functionOpen = False


            # Open a new function
            if lineTokens[0].lower() == "start":
                functionOpen = True
                functionCounter += 1
                functionLists.append([])

        return functionLists



#######################################################################################################
##  Create the total class modifiers from the full line of modifiers in the parseInputs directory structure
#######################################################################################################

    def getCompletedClass(self,className,functionContent):

        nextLine = "\n"
        indent = "    "
        classContent = ""
        classContent = nextLine+"class "+className+":"+nextLine
        classContent += indent+"\"class for modifying the "+className+" parameterList\""+nextLine

        classContent += nextLine+nextLine
        classContent += indent+"def __init__(self):"+nextLine
        classContent += indent+indent+"self.modifierName = \""+className+"\""+nextLine

        classContent += nextLine+nextLine
        classContent += indent+"def getName(self):"+nextLine
        classContent += indent+indent+"return self.modifierName"+nextLine

        classContent += nextLine+nextLine
        for fC in functionContent:
            classContent += indent+fC


        return classContent



#######################################################################################################
##  Create the modifier library
#######################################################################################################

    def createModifierLibrary(self):
        # First, create the parser class file
        nextLine = "\n"
        self.indent = "    "
        self.indent2 = self.indent+self.indent
        self.indent3 = self.indent2+self.indent
        self.indent4 = self.indent3+self.indent
        self.indent5 = self.indent4+self.indent
        self.indent6 = self.indent5+self.indent
        self.indent7 = self.indent6+self.indent
        self.indent8 = self.indent7+self.indent
        self.indent9 = self.indent8+self.indent
        self.indent10 = self.indent9+self.indent
        # First need to create init in the modifier directory
        initFile = open(self.modifierTargetPath+"/__init__.py","w+")
        initFile.write(" ")
        initFile.close()

        fileContent = nextLine
        for mLL in self.modifierLibraryList:
            fileContent += "from ."+mLL+" import *"+nextLine


        # class head
        fileContent += nextLine+nextLine
        fileContent += "class modifierLibrary:"+nextLine
        fileContent += self.indent+"\"modifier library\""+nextLine

        # init
        fileContent += nextLine+nextLine
        fileContent += self.indent+"def __init__(self):"+nextLine
        fileContent += self.indent2+"# create the modifier objects"+nextLine
        fileContent += self.indent2+"self.modifierList = []"+nextLine
        for mLL in self.modifierLibraryList:
            fileContent += self.indent2+"self.modifierList.append("+mLL+"())"+nextLine

        #execute modifiers
        fileContent += nextLine+nextLine
        fileContent += self.indent+"def executeModifiers(self,useModList,pLList):"+nextLine
        fileContent += self.indent2+"# Loop over the requested modifiers"+nextLine
        fileContent += self.indent2+"for uML in useModList:"+nextLine
        fileContent += self.indent3+"for mL in self.modifierList:"+nextLine
        fileContent += self.indent4+"if uML == mL.getName():"+nextLine
        fileContent += self.indent5+"pLList = mL.testForModification(pLList)"+nextLine

        # Return the modified parameter list
        fileContent += nextLine
        fileContent += self.indent2+"return pLList"+nextLine

        modifierLibrary = open(self.modifierTargetPath+"/charonModifierLibrary.py","w+")
        modifierLibrary.write(fileContent)
        modifierLibrary.close()


