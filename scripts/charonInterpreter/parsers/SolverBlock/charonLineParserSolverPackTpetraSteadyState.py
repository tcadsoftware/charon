
from __future__ import print_function
import copy


class charonLineParserSolverPackTpetraSteadyState:
    "SolverPackTpetraSteadyState parser"

    def __init__(self):
        # Register the parsing keys
        self.parserName = "SolverPackTpetraSteadyState"
        self.parsingKey = "use tpetra steady state solver pack"
        self.parsingKeyOptional = []
        self.interpreterHelpLine = "use tpetra steady state solver pack "
        self.interpreterQuickHelp = "Use solver settings generally optimal for use with Tpetra steady-state simulations."
        self.interpreterLongHelp = "Use solver settings generally optimal for use with Tpetra steady-state simulations."

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
        self.xmlDefaultLines.append("Charon->Solution Control,Compute Sensitivities,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control,Jacobian Operator,string,Have Jacobian")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction,Method,string,Newton")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton,Forcing Term Method,string,Constant")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton,Rescue Bad Newton Solve,bool,true")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Linear Solver,Tolerance,double,1.0e-6")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos,Linear Solver Type,string,Belos")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos,Solver Type,string,Block GMRES")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos->VerboseObject,Verbosity Level,string,none")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos->Solver Types->Block GMRES,Output Frequency,int,10")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos->Solver Types->Block GMRES,Output Style,int,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos->Solver Types->Block GMRES,Verbosity,int,33")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos->Solver Types->Block GMRES,Block Size,Maximum Iterations,int,300")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos->Solver Types->Block GMRES,Block Size,int,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos->Solver Types->Block GMRES,Num Blocks,int,20")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos->Solver Types->Block GMRES,Maximum Iterations,int,300")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Linear Solver Types->Belos->Solver Types->Block GMRES,Flexible Gmres,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos,Preconditioner Type,string,Ifpack2")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack2,Prec Type,string,SCHWARZ")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack2->VerboseObject,Verbosity Level,string,none")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack2->Ifpack2 Settings,schwarz: overlap level,int,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack2->Ifpack2 Settings,schwarz: use reordering,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Direction->Newton->Stratimikos Linear Solver->Stratimikos->Preconditioner Types->Ifpack2->Ifpack2 Settings,schwarz: combine mode,string,add")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX,Nonlinear Solver,string,Line Search Based")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Thyra Group Options,Function Scaling,string,Row Sum")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Line Search->Full Step,Full Step,double,1")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Line Search,Method,string,Full Step")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing,Output Precision,int,3")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing,Output Processor,int,0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Error,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Warning,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Outer Iteration,bool,true")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Parameters,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Details,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Linear Solver Details,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Stepper Iteration,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Stepper Details,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Printing->Output Information,Stepper Parameters,bool,false")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Solver Options,Status Test Check Type,string,Minimal")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests,Test Type,string,Combo")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests,Combo Type,string,OR")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests,Number of Tests,int,3")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0,Test Type,string,Combo")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0,Combo Type,string,AND")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0,Number of Tests,int,2")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0->Test 0,Test Type,string,NormF")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0->Test 0,Tolerance,double,1.0e-4")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0->Test 1,Test Type,string,NormWRMS")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0->Test 1,Tolerance,double,1.0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0->Test 1,Relative Tolerance,double,1.0e-09")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0->Test 1,Absolute Tolerance,double,1.0e-18")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0->Test 1,alpha,double,0.0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 0->Test 1,beta,double,1.0")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 1,Test Type,string,MaxIters")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 1,Maximum Iterations,int,30")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 2,Test Type,string,FiniteValue")
        self.xmlDefaultLines.append("Charon->Solution Control->NOX->Status Tests->Test 2,Vector Type,string,F Vector")
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
