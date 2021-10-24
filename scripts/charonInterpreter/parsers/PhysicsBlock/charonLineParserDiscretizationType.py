
from __future__ import print_function
import copy


class charonLineParserDiscretizationType:
    "DiscretizationType parser"

    def __init__(self):
        # Register the parsing keys
        self.parserName = "DiscretizationType"
        self.parsingKey = "standard discretization type"
        self.parsingKeyOptional = []
        self.parsingKeyOptional.append("drift diffusion gfem")
        self.parsingKeyOptional.append("drift diffusion effpg")
        self.parsingKeyOptional.append("nlp")
        self.parsingKeyOptional.append("drift diffusion cvfem")
        self.parsingKeyOptional.append("laplace gfem")
        self.parsingKeyOptional.append("laplace cvfem")
        self.parsingKeyOptional.append("ddlattice gfem")
        self.parsingKeyOptional.append("source term stabilization")
        self.parsingKeyOptional.append("lattice gfem")
        self.interpreterHelpLine = "standard discretization type is [ drift diffusion gfem [ drift diffusion effpg [ nlp [ drift diffusion cvfem [ laplace gfem [ laplace cvfem [ddlattice gfem [ with source term stabilization [ lattice gfem ]]]]]]]]] "
        self.interpreterQuickHelp = "Set the discretization type and the equations to be solved"
        self.interpreterLongHelp = "Set the discretization and the equations to be solved <> drift diffusion gfem sets the type to the stabilized Galerkin finite element method for the drift diffusion equations <> drift diffusion cvfem sets the type to the Scharfetter Gummel method for the drift diffusion equations <> drift diffusion effpg sets the type to the Exponentially-Fitted Flux Petrov-Glarkin stabilized Galerkin finite element method for the drift diffusion equations <> nlp sets the nonlinear poisson equations <> laplace gfem sets the type to the stabilized Galerkin finite element method for the Laplace equation <> laplace cvfem sets the type to a finite volume method for Laplace equation"

        # Register the xml required lines
        self.xmlRequiredLines = []
        self.xmlRequiredLinePriority = []
        self.xmlNewRequiredLines = []

        # Register the xml required arguments and their indexes
        self.xmlRequiredArgument = []
        self.xmlRequiredArgumentIndexes = []

        # Register the xml optional lines
        self.xmlOptionalLines = [[]]
        self.xmlOptionalLinePriority = [[]]
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Type,string,Drift Diffusion")
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Type,string,HGrad")
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Order,int,1")
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Integration Order,int,2")
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Solve Electron,string,True")
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Solve Hole,string,True")
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,SUPG Stabilization,string,On")
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Tau_E,string,Tanh")
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Tau_H,string,Tanh")
        self.xmlOptionalLines[0].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Length Scale,string,Stream")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[1].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Type,string,EFFPG Drift Diffusion")
        self.xmlOptionalLines[1].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Type,string,HGrad")
        self.xmlOptionalLines[1].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Order,int,1")
        self.xmlOptionalLines[1].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Integration Order,int,2")
        self.xmlOptionalLines[1].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Solve Electron,string,True")
        self.xmlOptionalLines[1].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Solve Hole,string,True")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[2].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Type,string,NLPoisson")
        self.xmlOptionalLines[2].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Type,string,HGrad")
        self.xmlOptionalLines[2].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Order,int,1")
        self.xmlOptionalLines[2].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Integration Order,int,2")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[3].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Type,string,SGCVFEM Drift Diffusion")
        self.xmlOptionalLines[3].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Type,string,HGrad")
        self.xmlOptionalLines[3].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Order,int,1")
        self.xmlOptionalLines[3].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Integration Order,int,2")
        self.xmlOptionalLines[3].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Solve Electron,string,True")
        self.xmlOptionalLines[3].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Solve Hole,string,True")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[4].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Type,string,Laplace")
        self.xmlOptionalLines[4].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Type,string,HGrad")
        self.xmlOptionalLines[4].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Order,int,1")
        self.xmlOptionalLines[4].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Integration Order,int,2")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[5].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Type,string,SGCVFEM Laplace")
        self.xmlOptionalLines[5].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Type,string,HGrad")
        self.xmlOptionalLines[5].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Order,int,1")
        self.xmlOptionalLines[5].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Integration Order,int,2")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Type,string,DDLattice")
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Type,string,HGrad")
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Order,int,1")
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Integration Order,int,2")
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Solve Electron,string,True")
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Solve Hole,string,True")
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,SUPG Stabilization,string,On")
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Tau_E,string,Tanh")
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Tau_H,string,Tanh")
        self.xmlOptionalLines[6].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Length Scale,string,Stream")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[7].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Add Source Term,string,On")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[8].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Type,string,Lattice")
        self.xmlOptionalLines[8].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Type,string,HGrad")
        self.xmlOptionalLines[8].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Basis Order,int,1")
        self.xmlOptionalLines[8].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS,Integration Order,int,2")
        self.xmlOptionalLines[8].append(" Charon->Physics Blocks->{physicsBlockName}->ANONYMOUS->Options,Heat Generation,string,Analytic")
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[2].append(2)
        self.xmlOptionalLinePriority[2].append(2)
        self.xmlOptionalLinePriority[2].append(2)
        self.xmlOptionalLinePriority[2].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[3].append(2)
        self.xmlOptionalLinePriority[3].append(2)
        self.xmlOptionalLinePriority[3].append(2)
        self.xmlOptionalLinePriority[3].append(2)
        self.xmlOptionalLinePriority[3].append(2)
        self.xmlOptionalLinePriority[3].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[4].append(2)
        self.xmlOptionalLinePriority[4].append(2)
        self.xmlOptionalLinePriority[4].append(2)
        self.xmlOptionalLinePriority[4].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[5].append(2)
        self.xmlOptionalLinePriority[5].append(2)
        self.xmlOptionalLinePriority[5].append(2)
        self.xmlOptionalLinePriority[5].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority[6].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[7].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[8].append(2)
        self.xmlOptionalLinePriority[8].append(2)
        self.xmlOptionalLinePriority[8].append(2)
        self.xmlOptionalLinePriority[8].append(2)
        self.xmlOptionalLinePriority[8].append(2)

        # Register the xml optional arguments and their indexes
        self.xmlOptionalArgument = [[], [], [], [], [], [], [], [], []]
        self.xmlOptionalArgumentIndexes = [[], [], [], [], [], [], [], [], []]

        # Register the xml default lines
        self.xmlDefaultLines = []
        self.xmlDefaultLinePriority = []

        self.xmlReturned = []
        self.xmlPriorityCode = []



    def isThisMe(self,tokenizer,line):
        # Tokenize the line
        lineTokens = tokenizer.tokenize(line)
        # Tokenize the parsing key
        parsingTokens = self.parsingKey.split()
        returnType = True
        for itoken in range(len(parsingTokens)):
            if itoken+1 > len(lineTokens):
                return False
            if lineTokens[itoken].lower() != parsingTokens[itoken].lower():
                returnType = False
        return returnType



    def getName(self):
        # Return parser name
         return self.parserName



    def getHelp(self,verbosity):
        # Return help content
        if verbosity.lower() == "long":
            return (self.interpreterHelpLine,self.interpreterLongHelp)
        else:
            return (self.interpreterHelpLine,self.interpreterQuickHelp)



    def generateXML(self,tokenizer,line):
        # Tokenize the line
        lineTokens = tokenizer.tokenize(line)
        self.xmlNewRequiredLines[:] = []
        for xL in self.xmlRequiredLines:
            self.xmlNewRequiredLines.append(xL)
        for ipar in range(len(self.xmlRequiredArgument)):
            line.replace(self.xmlRequiredArgument[ipar],lineTokens[int(self.xmlRequiredArgumentIndexes[ipar])])
            for iRLine in range(len(self.xmlRequiredLines)):
                self.xmlNewRequiredLines[iRLine]=self.xmlNewRequiredLines[iRLine].replace(self.xmlRequiredArgument[ipar],lineTokens[int(self.xmlRequiredArgumentIndexes[ipar])])
        for index,xmlLine in enumerate(self.xmlNewRequiredLines):
            self.xmlReturned.append(xmlLine)
            self.xmlPriorityCode.append(self.xmlRequiredLinePriority[index]) #required lines have priority code 2
        # Look over input line to see if any options are called out.
        optCounter = 0
        optIndex = 0
        for optKey in self.parsingKeyOptional:
            # Tokenize the opt keys
            foundOptionalKey = False
            optKeyTokens = optKey.split()
            for iLT in range(len(lineTokens)):
                if lineTokens[iLT].lower() == optKeyTokens[0]:
                    if len(optKeyTokens) == 1:
                        optIndex = iLT
                        foundOptionalKey = True
                    else:
                        for iPK in range(len(optKeyTokens)-1):
                            optIndex = iLT
                            if iLT+iPK+1 > len(lineTokens)-1:
                                continue
                            if optKeyTokens[iPK+1] == lineTokens[iLT+iPK+1].lower():
                                if iPK+2 == len(optKeyTokens):
                                    foundOptionalKey = True
                                else:
                                    continue
            #Found the key, now create the xml line
            if foundOptionalKey == True:
                self.Returned=copy.deepcopy(self.xmlOptionalLines[optCounter])
                for iopt in range(len(self.xmlOptionalLines[optCounter])):
                    for ipar in range(len(self.xmlOptionalArgument[optCounter])):
                        self.Returned[iopt] = self.Returned[iopt].replace(self.xmlOptionalArgument[optCounter][ipar],lineTokens[optIndex+int(self.xmlOptionalArgumentIndexes[optCounter][ipar])])
                    for ipar in range(len(self.xmlRequiredArgument)):
                        self.Returned[iopt] = self.Returned[iopt].replace(self.xmlRequiredArgument[ipar],lineTokens[int(self.xmlRequiredArgumentIndexes[ipar])])
                    self.xmlReturned.append(self.Returned[iopt])
                    self.xmlPriorityCode.append(2) #optional lines have priority code 2
            optCounter += 1
        for xmlLine in self.xmlDefaultLines:
            self.xmlReturned.append(xmlLine)
            self.xmlPriorityCode.append(1) #optional lines have priority code 1

        return (self.xmlReturned,self.xmlPriorityCode)
