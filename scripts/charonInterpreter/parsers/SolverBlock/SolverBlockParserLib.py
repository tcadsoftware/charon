
try:
    import coloramaDISABLED as colors
except ImportError:
    class stubColors:
        "subs for colors when colors doesn't exist on system"
    
        def __init__(self):
            self.Fore = colorClass()
            self.Back = colorClass()
            self.Style = styleClass()
    
    class colorClass():
        "stubbed color class"
    
        def __init__(self):
            self.BLACK = ""
            self.BLUE = ""
            self.WHITE = ""
            self.RED = ""
            self.GREEN = ""
    
    class styleClass():
        "stubbed style class"
    
        def __init__(self):
            self.RESET_ALL = ""
    colors = stubColors()

import sys
from .charonLineParserSolverPack5 import *
from .charonLineParserSolverPack3 import *
from .charonLineParserSolverPack1 import *
from .charonLineParserSolverPack4 import *
from .charonLineParserLineSearchMethod import *
from .charonLineParserThyraFunctionScaling import *
from .charonLineParserTransientSolverPack2 import *
from .charonLineParserIfpackLevelOfFill import *
from .charonLineParserTimeIntegratorMaxTimeStepSize import *
from .charonLineParserTimeIntegratorFinalTime import *
from .charonLineParserLinearSolverType import *
from .charonLineParserTimeIntegratorFixedStepSize import *
from .charonLineParserMaxNonlinearIterations import *
from .charonLineParserTransientSolverPack1 import *
from .charonLineParserIfpackReorderingType import *
from .charonLineParserSolverPack2 import *
from .charonLineParserTimeIntegratorMinTimeStepSize import *
from .charonLineParserSolverPackHB import *
from .charonLineParserAztecMaxIterations import *
from .charonLineParserSolverPack6 import *
from .charonLineParserGMRESILU import *
from .charonLineParserTimeIntegratorInitialStepSize import *
from .charonLineParserTimeIntegratorVariableStepSize import *
from .charonLineParserSolverPack7 import *
from .charonLineParserTimeIntegratorAbsoluteTolerance import *
from .charonLineParserTimeIntegratorRelativeTolerance import *
from .charonLineParserSolverTolerance import *
from .charonLineParserNoxSolver import *
from .charonLineParserLineSearchFullStep import *
from .charonLineParserAztecTolerance import *
from .charonLineParserAztecKrylovSubSize import *
from .charonBlockParserWRMSParameters import *
from .WRMSParameters.WRMSParametersParserLib import *



class SolverBlockParserLib:
    "This is the  SolverBlockParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "SolverBlockParserLib"
        # create the linparser objects 
        self.lineParsers = []
        self.lineParsers.append(charonLineParserSolverPack5())
        self.lineParsers.append(charonLineParserSolverPack3())
        self.lineParsers.append(charonLineParserSolverPack1())
        self.lineParsers.append(charonLineParserSolverPack4())
        self.lineParsers.append(charonLineParserLineSearchMethod())
        self.lineParsers.append(charonLineParserThyraFunctionScaling())
        self.lineParsers.append(charonLineParserTransientSolverPack2())
        self.lineParsers.append(charonLineParserIfpackLevelOfFill())
        self.lineParsers.append(charonLineParserTimeIntegratorMaxTimeStepSize())
        self.lineParsers.append(charonLineParserTimeIntegratorFinalTime())
        self.lineParsers.append(charonLineParserLinearSolverType())
        self.lineParsers.append(charonLineParserTimeIntegratorFixedStepSize())
        self.lineParsers.append(charonLineParserMaxNonlinearIterations())
        self.lineParsers.append(charonLineParserTransientSolverPack1())
        self.lineParsers.append(charonLineParserIfpackReorderingType())
        self.lineParsers.append(charonLineParserSolverPack2())
        self.lineParsers.append(charonLineParserTimeIntegratorMinTimeStepSize())
        self.lineParsers.append(charonLineParserSolverPackHB())
        self.lineParsers.append(charonLineParserAztecMaxIterations())
        self.lineParsers.append(charonLineParserSolverPack6())
        self.lineParsers.append(charonLineParserGMRESILU())
        self.lineParsers.append(charonLineParserTimeIntegratorInitialStepSize())
        self.lineParsers.append(charonLineParserTimeIntegratorVariableStepSize())
        self.lineParsers.append(charonLineParserSolverPack7())
        self.lineParsers.append(charonLineParserTimeIntegratorAbsoluteTolerance())
        self.lineParsers.append(charonLineParserTimeIntegratorRelativeTolerance())
        self.lineParsers.append(charonLineParserSolverTolerance())
        self.lineParsers.append(charonLineParserNoxSolver())
        self.lineParsers.append(charonLineParserLineSearchFullStep())
        self.lineParsers.append(charonLineParserAztecTolerance())
        self.lineParsers.append(charonLineParserAztecKrylovSubSize())
        # create the blockparser objects 
        self.blockParsers = []
        self.blockParsers.append([charonBlockParserWRMSParameters(),WRMSParametersParserLib()])
        # create the parserLibrary objects 
        parserLibraries = []
        parserLibraries.append(WRMSParametersParserLib())


    def isThisMyLine(self,tokenizer,line):
        for lP in self.lineParsers:
            self.isThisMe = lP.isThisMe(tokenizer,line)
            if self.isThisMe == True:
                return (True,lP)
        return (False,None)


    def isThisMyBlock(self,tokenizer,line):
        for bP in self.blockParsers:
            self.isThisMe = bP[0].isThisMe(tokenizer,line)
            if self.isThisMe == True:
                return (True,bP[0],bP[1])
        return (False,None,None)


    def generateHelp(self,genHelp,indent):
        self.addIndent = "     "
        cRStyle = ""
        for lP in self.lineParsers:
            (self.helpLine,self.helpContent) = lP.getHelp(genHelp)
            self.helpContentList = self.helpContent.split("<>")
            print (cRStyle+indent+colors.Fore.RED+colors.Back.WHITE+self.helpLine)
            cRStyle = "\n"
            for hCL in self.helpContentList:
                print ("\t"+indent+colors.Fore.BLUE+colors.Back.WHITE+hCL.lstrip())
        for bP in range(len(self.blockParsers)):
            print (indent+colors.Fore.GREEN+colors.Back.WHITE+self.blockParsers[bP][0].getHelpLine().lstrip())
            self.blockParsers[bP][1].generateHelp(genHelp,indent+self.addIndent)
            print (indent+colors.Fore.GREEN+colors.Back.WHITE+self.blockParsers[bP][0].getHelpLine().replace("start","end").lstrip())
            print (indent+colors.Style.RESET_ALL)


    def getName(self):
        return self.parserLibName
