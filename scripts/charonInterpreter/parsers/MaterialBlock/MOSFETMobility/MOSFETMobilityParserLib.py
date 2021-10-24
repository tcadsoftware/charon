
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
from .charonBlockParserHoleShirahataMobility import *
from .charonBlockParserElectronShirahataMobility import *
from .charonBlockParserHoleBulkPhilipsThomasMobility import *
from .charonBlockParserElectronBulkPhilipsThomasMobility import *
from .HoleShirahataMobility.HoleShirahataMobilityParserLib import *
from .ElectronShirahataMobility.ElectronShirahataMobilityParserLib import *
from .HoleBulkPhilipsThomasMobility.HoleBulkPhilipsThomasMobilityParserLib import *
from .ElectronBulkPhilipsThomasMobility.ElectronBulkPhilipsThomasMobilityParserLib import *



class MOSFETMobilityParserLib:
    "This is the  MOSFETMobilityParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "MOSFETMobilityParserLib"
        # create the linparser objects 
        self.lineParsers = []
        # create the blockparser objects 
        self.blockParsers = []
        self.blockParsers.append([charonBlockParserHoleShirahataMobility(),HoleShirahataMobilityParserLib()])
        self.blockParsers.append([charonBlockParserElectronShirahataMobility(),ElectronShirahataMobilityParserLib()])
        self.blockParsers.append([charonBlockParserHoleBulkPhilipsThomasMobility(),HoleBulkPhilipsThomasMobilityParserLib()])
        self.blockParsers.append([charonBlockParserElectronBulkPhilipsThomasMobility(),ElectronBulkPhilipsThomasMobilityParserLib()])
        # create the parserLibrary objects 
        parserLibraries = []
        parserLibraries.append(HoleShirahataMobilityParserLib())
        parserLibraries.append(ElectronShirahataMobilityParserLib())
        parserLibraries.append(HoleBulkPhilipsThomasMobilityParserLib())
        parserLibraries.append(ElectronBulkPhilipsThomasMobilityParserLib())


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
