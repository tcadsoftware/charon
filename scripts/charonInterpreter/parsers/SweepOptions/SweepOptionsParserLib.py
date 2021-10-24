
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
from .charonLineParserContinuationMethod import *
from .charonLineParserSuccessfulStepIncreaseFac import *
from .charonLineParserMinValue import *
from .charonLineParserMinimumStepSize import *
from .charonLineParserMaxValue import *
from .charonLineParserInitialStepSize import *
from .charonLineParserFailedStepReductionFac import *
from .charonLineParserMaximumStepSize import *
from .charonLineParserPredictorMethod import *
from .charonLineParserMaximumNumberSteps import *
from .charonLineParserStepSizeAggressiveness import *



class SweepOptionsParserLib:
    "This is the  SweepOptionsParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "SweepOptionsParserLib"
        # create the linparser objects 
        self.lineParsers = []
        self.lineParsers.append(charonLineParserContinuationMethod())
        self.lineParsers.append(charonLineParserSuccessfulStepIncreaseFac())
        self.lineParsers.append(charonLineParserMinValue())
        self.lineParsers.append(charonLineParserMinimumStepSize())
        self.lineParsers.append(charonLineParserMaxValue())
        self.lineParsers.append(charonLineParserInitialStepSize())
        self.lineParsers.append(charonLineParserFailedStepReductionFac())
        self.lineParsers.append(charonLineParserMaximumStepSize())
        self.lineParsers.append(charonLineParserPredictorMethod())
        self.lineParsers.append(charonLineParserMaximumNumberSteps())
        self.lineParsers.append(charonLineParserStepSizeAggressiveness())
        # create the blockparser objects 
        self.blockParsers = []
        # create the parserLibrary objects 
        parserLibraries = []


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
