
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
from .charonLineParserRadiativeRecombination import *
from .charonLineParserIncompleteIonizationAcceptor import *
from .charonLineParserSRHRecombination import *
from .charonLineParserAugerRecombination import *
from .charonLineParserTrapSRH import *
from .charonLineParserPecletLengthScale import *
from .charonLineParserSUPGStabilizationToggle import *
from .charonLineParserDefectCluster import *
from .charonLineParserDiscretizationMethod import *
from .charonLineParserSetMaterialModel import *
from .charonLineParserIncompleteIonizationDonor import *
from .charonLineParserEmpiricalModel import *
from .charonLineParserDiscretizationType import *
from .charonLineParserFermiDirac import *
from .charonLineParserBulkFixedCharge import *
from .charonLineParserAvalancheGen import *
from .charonLineParserDrivingForce import *
from .charonLineParserBandGapNarrowing import *
from .charonLineParserParticleStrike import *
from .charonLineParserSetGeometryBlock import *
from .charonLineParserHeatGeneration import *
from .charonLineParserTrapCharge import *
from .charonLineParserDiscontinuity import *
from .charonLineParserHarmonicBalanceAnalysis import *
from .charonLineParserOpticalGeneration import *
from .charonBlockParserHarmonicBalanceBlock import *
from .HarmonicBalanceBlock.HarmonicBalanceBlockParserLib import *



class PhysicsBlockParserLib:
    "This is the  PhysicsBlockParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "PhysicsBlockParserLib"
        # create the linparser objects 
        self.lineParsers = []
        self.lineParsers.append(charonLineParserRadiativeRecombination())
        self.lineParsers.append(charonLineParserIncompleteIonizationAcceptor())
        self.lineParsers.append(charonLineParserSRHRecombination())
        self.lineParsers.append(charonLineParserAugerRecombination())
        self.lineParsers.append(charonLineParserTrapSRH())
        self.lineParsers.append(charonLineParserPecletLengthScale())
        self.lineParsers.append(charonLineParserSUPGStabilizationToggle())
        self.lineParsers.append(charonLineParserDefectCluster())
        self.lineParsers.append(charonLineParserDiscretizationMethod())
        self.lineParsers.append(charonLineParserSetMaterialModel())
        self.lineParsers.append(charonLineParserIncompleteIonizationDonor())
        self.lineParsers.append(charonLineParserEmpiricalModel())
        self.lineParsers.append(charonLineParserDiscretizationType())
        self.lineParsers.append(charonLineParserFermiDirac())
        self.lineParsers.append(charonLineParserBulkFixedCharge())
        self.lineParsers.append(charonLineParserAvalancheGen())
        self.lineParsers.append(charonLineParserDrivingForce())
        self.lineParsers.append(charonLineParserBandGapNarrowing())
        self.lineParsers.append(charonLineParserParticleStrike())
        self.lineParsers.append(charonLineParserSetGeometryBlock())
        self.lineParsers.append(charonLineParserHeatGeneration())
        self.lineParsers.append(charonLineParserTrapCharge())
        self.lineParsers.append(charonLineParserDiscontinuity())
        self.lineParsers.append(charonLineParserHarmonicBalanceAnalysis())
        self.lineParsers.append(charonLineParserOpticalGeneration())
        # create the blockparser objects 
        self.blockParsers = []
        self.blockParsers.append([charonBlockParserHarmonicBalanceBlock(),HarmonicBalanceBlockParserLib()])
        # create the parserLibrary objects 
        parserLibraries = []
        parserLibraries.append(HarmonicBalanceBlockParserLib())


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
