
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
from .charonLineParserNetlistTemplate import *
from .charonLineParserCascadeDensity import *
from .charonLineParserShepardPower import *
from .charonLineParserInputFileType import *
from .charonLineParserSettleTime import *
from .charonLineParserClusterLocations import *
from .charonLineParserOnsetTime import *
from .charonLineParserNumberOfInputFiles import *
from .charonLineParserPulseFile import *
from .charonLineParserVerbosityLevel import *
from .charonLineParserSingleFileClusters import *
from .charonLineParserInterpolantMethod import *
from .charonLineParserSaveClusterCoefficients import *



class DefectClusterParametersParserLib:
    "This is the  DefectClusterParametersParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "DefectClusterParametersParserLib"
        # create the linparser objects 
        self.lineParsers = []
        self.lineParsers.append(charonLineParserNetlistTemplate())
        self.lineParsers.append(charonLineParserCascadeDensity())
        self.lineParsers.append(charonLineParserShepardPower())
        self.lineParsers.append(charonLineParserInputFileType())
        self.lineParsers.append(charonLineParserSettleTime())
        self.lineParsers.append(charonLineParserClusterLocations())
        self.lineParsers.append(charonLineParserOnsetTime())
        self.lineParsers.append(charonLineParserNumberOfInputFiles())
        self.lineParsers.append(charonLineParserPulseFile())
        self.lineParsers.append(charonLineParserVerbosityLevel())
        self.lineParsers.append(charonLineParserSingleFileClusters())
        self.lineParsers.append(charonLineParserInterpolantMethod())
        self.lineParsers.append(charonLineParserSaveClusterCoefficients())
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
