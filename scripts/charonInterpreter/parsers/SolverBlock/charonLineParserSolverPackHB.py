
from __future__ import print_function
import copy


class charonLineParserSolverPackHB:
    "SolverPackHB parser"

    def __init__(self):
        # Register the parsing keys
        self.parserName = "SolverPackHB"
        self.parsingKey = "use solver pack hb"
        self.parsingKeyOptional = []
        self.interpreterHelpLine = "use solver pack HB "
        self.interpreterQuickHelp = "Specify the use of the frequency domain solver pack HB"
        self.interpreterLongHelp = "Specify the use of the frequency domain solver pack HB"

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
        self.xmlDefaultLines.append("Charon->Solution Control,Piro Solver,string,NOX")
        self.xmlDefaultLines.append("Charon->Solution Control,Compute Sensitivities,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control,Jacobian Operator,string,Have Jacobian")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX,Nonlinear Solver,string,Line Search Based")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction,Method,string,Newton")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton,Forcing Term Method,string,Constant")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton,Rescue Bad Newton Solve,bool,true")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos,Linear Solver Type,string,AztecOO")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos,Preconditioner Type,string,Ifpack")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve,Max Iterations,int,300")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve,Tolerance,double,1e-5")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve->AztecOO Settings,Aztec Solver,string,GMRES")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve->AztecOO Settings,Convergence Test,string,r0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve->AztecOO Settings,Size of Krylov Subspace,int,700")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->AztecOO->Forward Solve->AztecOO Settings,Output Frequency,int,10")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack,Overlap,int,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack,Prec Type,string,point relaxation")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack->Ifpack Settings,relaxation: type,string,Jacobi")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack->Ifpack Settings,fact: level-of-fill,int,10")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack->Ifpack Settings,schwarz: reordering type,string,rcm")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Line Search,Method,string,Polynomial")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Line Search->Polynomial,Maximum Iterations,int,10")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Line Search->Polynomial,Recovery Step Type,string,Last Computed Step")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing,Output Precision,int,3")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing,Output Processor,int,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Error,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Warning,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Outer Iteration,bool,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Parameters,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Details,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Linear Solver Details,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Stepper Iteration,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Stepper Details,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Stepper Parameters,bool,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Solver Options,Status Test Check Type,string,Minimal")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests,Test Type,string,Combo")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests,Combo Type,string,OR")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests,Number of Tests,int,2")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0,Test Type,string,NormF")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0,Tolerance,double,1.0e-6")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 1,Test Type,string,MaxIters")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 1,Maximum Iterations,int,30")
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
