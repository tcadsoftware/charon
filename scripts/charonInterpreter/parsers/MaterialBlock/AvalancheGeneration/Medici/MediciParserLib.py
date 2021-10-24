
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
from .charonLineParserHole_a2 import *
from .charonLineParserElectron_a0 import *
from .charonLineParserHole_hbaromega import *
from .charonLineParserElectron_delta import *
from .charonLineParserHole_lambda300 import *
from .charonLineParserElectron_lambda300 import *
from .charonLineParserElectron_hbaromega import *
from .charonLineParserHole_a1 import *
from .charonLineParserElectron_a1 import *
from .charonLineParserElectron_E0 import *
from .charonLineParserHole_delta import *
from .charonLineParserCriticalField import *
from .charonLineParserElectron_a2 import *
from .charonLineParserHole_a0 import *
from .charonLineParserHole_E0 import *



class MediciParserLib:
    "This is the  MediciParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "MediciParserLib"
        # create the linparser objects 
        self.lineParsers = []
        self.lineParsers.append(charonLineParserHole_a2())
        self.lineParsers.append(charonLineParserElectron_a0())
        self.lineParsers.append(charonLineParserHole_hbaromega())
        self.lineParsers.append(charonLineParserElectron_delta())
        self.lineParsers.append(charonLineParserHole_lambda300())
        self.lineParsers.append(charonLineParserElectron_lambda300())
        self.lineParsers.append(charonLineParserElectron_hbaromega())
        self.lineParsers.append(charonLineParserHole_a1())
        self.lineParsers.append(charonLineParserElectron_a1())
        self.lineParsers.append(charonLineParserElectron_E0())
        self.lineParsers.append(charonLineParserHole_delta())
        self.lineParsers.append(charonLineParserCriticalField())
        self.lineParsers.append(charonLineParserElectron_a2())
        self.lineParsers.append(charonLineParserHole_a0())
        self.lineParsers.append(charonLineParserHole_E0())
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
