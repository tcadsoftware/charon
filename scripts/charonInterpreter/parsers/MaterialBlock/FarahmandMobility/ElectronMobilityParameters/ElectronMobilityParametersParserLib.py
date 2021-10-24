
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
from .charonLineParserHighFieldMobility import *
from .charonLineParserEps import *
from .charonLineParserAn import *
from .charonLineParserMu1 import *
from .charonLineParserNcrit import *
from .charonLineParserGamma import *
from .charonLineParserVsat import *
from .charonLineParserDrivingForce import *
from .charonLineParserN1 import *
from .charonLineParserDelta import *
from .charonLineParserBeta import *
from .charonLineParserEc import *
from .charonLineParserAlpha import *
from .charonLineParserMu2 import *
from .charonLineParserN2 import *



class ElectronMobilityParametersParserLib:
    "This is the  ElectronMobilityParametersParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "ElectronMobilityParametersParserLib"
        # create the linparser objects 
        self.lineParsers = []
        self.lineParsers.append(charonLineParserHighFieldMobility())
        self.lineParsers.append(charonLineParserEps())
        self.lineParsers.append(charonLineParserAn())
        self.lineParsers.append(charonLineParserMu1())
        self.lineParsers.append(charonLineParserNcrit())
        self.lineParsers.append(charonLineParserGamma())
        self.lineParsers.append(charonLineParserVsat())
        self.lineParsers.append(charonLineParserDrivingForce())
        self.lineParsers.append(charonLineParserN1())
        self.lineParsers.append(charonLineParserDelta())
        self.lineParsers.append(charonLineParserBeta())
        self.lineParsers.append(charonLineParserEc())
        self.lineParsers.append(charonLineParserAlpha())
        self.lineParsers.append(charonLineParserMu2())
        self.lineParsers.append(charonLineParserN2())
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
