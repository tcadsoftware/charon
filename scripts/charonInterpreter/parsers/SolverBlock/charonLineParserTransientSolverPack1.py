
from __future__ import print_function
import copy


class charonLineParserTransientSolverPack1:
    "TransientSolverPack1 parser"

    def __init__(self):
        # Register the parsing keys
        self.parserName = "TransientSolverPack1"
        self.parsingKey = "use transient solver pack 1"
        self.parsingKeyOptional = []
        self.interpreterHelpLine = "use transient solver pack 1 "
        self.interpreterQuickHelp = "Specify the use of the transient solver pack 1"
        self.interpreterLongHelp = "Specify the use of the transient solver pack 1"

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

        # Register the xml optional arguments and their indexes
        self.xmlOptionalArgument = []
        self.xmlOptionalArgumentIndexes = []

        # Register the xml default lines
        self.xmlDefaultLines = []
        self.xmlDefaultLinePriority = []
        self.xmlDefaultLines.append("Charon->Solution Control,Piro Solver,string,Rythmos")
        self.xmlDefaultLines.append("Charon->Solution Control,Compute Sensitivities,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control,Jacobian Operator,string,Have Jacobian")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction,Method,string,Newton")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton,Forcing Term Method,string,Constant")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton,Rescue Bad Newton Solve,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos,Linear Solver Type,string,AztecOO")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve->AztecOO Settings,Aztec Solver,string,GMRES")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve->AztecOO Settings,Convergence Test,string,r0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve->AztecOO Settings,Size of Krylov Subspace,int,200")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve->AztecOO Settings,Output Frequency,int,10")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve,Max Iterations,int,300")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve,Tolerance,double,1e-5")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos,Preconditioner Type,string,Ifpack")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack,Overlap,int,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack,Prec Type,string,ILU")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack->Ifpack Settings,fact: drop tolerance,double,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack->Ifpack Settings,fact: ilut level-of-fill,double,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack->Ifpack Settings,fact: level-of-fill,int,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Line Search->Full Step,Full Step,double,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Line Search,Method,string,Backtrack")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX,Nonlinear Solver,string,Line Search Based")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing,Output Precision,int,3")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing,Output Processor,int,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Error,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Warning,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Outer Iteration,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Parameters,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Details,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Linear Solver Details,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Stepper Iteration,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Stepper Details,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Stepper Parameters,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Solver Options,Status Test Check Type,string,Minimal")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests,Test Type,string,Combo")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests,Combo Type,string,OR")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests,Number of Tests,int,2")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0,Test Type,string,NormF")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0,Tolerance,double,1.0e-8")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 1,Test Type,string,MaxIters")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 1,Maximum Iterations,int,10")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos,Nonlinear Solver Type,string,NOX")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos,Final Time,double,5.5e-4")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos,Stepper Type,string,BDF")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos,Step Control Strategy Type,string,ImplicitBDFRamping")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos,Rythmos Integration Control Strategy,string,Simple")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Stepper->VerboseObject,Verbosity Level,string,none")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Integrator->VerboseObject,Verbosity Level,string,none")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Integration Control,Take Variable Steps,bool,true")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Number of Constant First Order Steps,int,5")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Initial Step Size,double,1e-5")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Min Step Size,double,1e-7")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Max Step Size,double,1.0e-4")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Step Size Decrease Factor,double,0.2")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Step Size Increase Factor,double,1.5")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Min Order,int,1")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Max Order,int,1")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Absolute Error Tolerance,double,1.0e-6")
        self.xmlDefaultLines.append("Charon->Solution Control->Rythmos->Rythmos Step Control Strategy,Relative Error Tolerance,double,1.0e-3")
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)
        self.xmlDefaultLinePriority.append(1)

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
