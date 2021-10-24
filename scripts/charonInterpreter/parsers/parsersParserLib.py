
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
from .charonLineParserContactOnInsulator import *
from .charonLineParserConstantBC import *
from .charonLineParserBlockPhysicsMap import *
from .charonLineParserWriteLinearSystem import *
from .charonLineParserFixedLatticeTemperature import *
from .charonLineParserThermalBC import *
from .charonLineParserAssemblyWorkset import *
from .charonLineParser1DBJTBase import *
from .charonLineParserConcentrationScaling import *
from .charonLineParserImportStateFile import *
from .charonLineParserInputMeshScaleFactor import *
from .charonLineParserDirectXMLInput import *
from .charonLineParserTimePeriodicBC import *
from .charonLineParserNeumannThermalBC import *
from .charonLineParserLinearRampBC import *
from .charonLineParserResistorContact import *
from .charonLineParserInitialConditions import *
from .charonLineParserOutputTimings import *
from .charonLineParserTestLine import *
from .charonLineParserOhmicBC import *
from .charonLineParserInitialConditionsHarmonicBalance import *
from .charonLineParserScaleInitialConditions import *
from .charonLineParserMmsBC import *
from .charonLineParserConstantCurrentContact import *
from .charonBlockParserEmpiricalModelParameters import *
from .charonBlockParserPhysicsBlock import *
from .charonBlockParserHarmonicBalanceBC import *
from .charonBlockParserSolverBlock import *
from .charonBlockParserHeterojunction import *
from .charonBlockParserMaterialBlock import *
from .charonBlockParserSurfaceChargeBC import *
from .charonBlockParserSweepOptions import *
from .charonBlockParserSchottkyBC import *
from .charonBlockParserOutputParameters import *
from .EmpiricalModelParameters.EmpiricalModelParametersParserLib import *
from .PhysicsBlock.PhysicsBlockParserLib import *
from .HarmonicBalanceBC.HarmonicBalanceBCParserLib import *
from .SolverBlock.SolverBlockParserLib import *
from .Heterojunction.HeterojunctionParserLib import *
from .MaterialBlock.MaterialBlockParserLib import *
from .SurfaceChargeBC.SurfaceChargeBCParserLib import *
from .SweepOptions.SweepOptionsParserLib import *
from .SchottkyBC.SchottkyBCParserLib import *
from .OutputParameters.OutputParametersParserLib import *



class parsersParserLib:
    "This is the  parsersParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "parsersParserLib"
        # create the linparser objects 
        self.lineParsers = []
        self.lineParsers.append(charonLineParserContactOnInsulator())
        self.lineParsers.append(charonLineParserConstantBC())
        self.lineParsers.append(charonLineParserBlockPhysicsMap())
        self.lineParsers.append(charonLineParserWriteLinearSystem())
        self.lineParsers.append(charonLineParserFixedLatticeTemperature())
        self.lineParsers.append(charonLineParserThermalBC())
        self.lineParsers.append(charonLineParserAssemblyWorkset())
        self.lineParsers.append(charonLineParser1DBJTBase())
        self.lineParsers.append(charonLineParserConcentrationScaling())
        self.lineParsers.append(charonLineParserImportStateFile())
        self.lineParsers.append(charonLineParserInputMeshScaleFactor())
        self.lineParsers.append(charonLineParserDirectXMLInput())
        self.lineParsers.append(charonLineParserTimePeriodicBC())
        self.lineParsers.append(charonLineParserNeumannThermalBC())
        self.lineParsers.append(charonLineParserLinearRampBC())
        self.lineParsers.append(charonLineParserResistorContact())
        self.lineParsers.append(charonLineParserInitialConditions())
        self.lineParsers.append(charonLineParserOutputTimings())
        self.lineParsers.append(charonLineParserTestLine())
        self.lineParsers.append(charonLineParserOhmicBC())
        self.lineParsers.append(charonLineParserInitialConditionsHarmonicBalance())
        self.lineParsers.append(charonLineParserScaleInitialConditions())
        self.lineParsers.append(charonLineParserMmsBC())
        self.lineParsers.append(charonLineParserConstantCurrentContact())
        # create the blockparser objects 
        self.blockParsers = []
        self.blockParsers.append([charonBlockParserEmpiricalModelParameters(),EmpiricalModelParametersParserLib()])
        self.blockParsers.append([charonBlockParserPhysicsBlock(),PhysicsBlockParserLib()])
        self.blockParsers.append([charonBlockParserHarmonicBalanceBC(),HarmonicBalanceBCParserLib()])
        self.blockParsers.append([charonBlockParserSolverBlock(),SolverBlockParserLib()])
        self.blockParsers.append([charonBlockParserHeterojunction(),HeterojunctionParserLib()])
        self.blockParsers.append([charonBlockParserMaterialBlock(),MaterialBlockParserLib()])
        self.blockParsers.append([charonBlockParserSurfaceChargeBC(),SurfaceChargeBCParserLib()])
        self.blockParsers.append([charonBlockParserSweepOptions(),SweepOptionsParserLib()])
        self.blockParsers.append([charonBlockParserSchottkyBC(),SchottkyBCParserLib()])
        self.blockParsers.append([charonBlockParserOutputParameters(),OutputParametersParserLib()])
        # create the parserLibrary objects 
        parserLibraries = []
        parserLibraries.append(EmpiricalModelParametersParserLib())
        parserLibraries.append(PhysicsBlockParserLib())
        parserLibraries.append(HarmonicBalanceBCParserLib())
        parserLibraries.append(SolverBlockParserLib())
        parserLibraries.append(HeterojunctionParserLib())
        parserLibraries.append(MaterialBlockParserLib())
        parserLibraries.append(SurfaceChargeBCParserLib())
        parserLibraries.append(SweepOptionsParserLib())
        parserLibraries.append(SchottkyBCParserLib())
        parserLibraries.append(OutputParametersParserLib())


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
