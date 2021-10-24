
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
from .charonLineParserFixedMobility import *
from .charonLineParserAffinity import *
from .charonLineParserMaterialName import *
from .charonLineParserRadiativeRecombinationCoefficient import *
from .charonLineParserDopingFile import *
from .charonLineParserReferenceMaterial import *
from .charonLineParserIntrinsicConcentrationModel import *
from .charonLineParserRelativePermittivity import *
from .charonLineParserBandgap import *
from .charonLineParserMMSDoping import *
from .charonBlockParserAroraMobility import *
from .charonBlockParserGlobalMMSParameters import *
from .charonBlockParserStepDoping import *
from .charonBlockParserMasettiMobility import *
from .charonBlockParserLucentMobility import *
from .charonBlockParserHeatGeneration import *
from .charonBlockParserCarrierLifetime import *
from .charonBlockParserMGaussDoping import *
from .charonBlockParserBandGapNarrowing import *
from .charonBlockParserIncompleteIonizationDonor import *
from .charonBlockParserTrapSRH import *
from .charonBlockParserIncompleteIonizationAcceptor import *
from .charonBlockParserBulkFixedChargeParameters import *
from .charonBlockParserOpticalGeneration import *
from .charonBlockParserEffectiveDOS import *
from .charonBlockParserUniformDoping import *
from .charonBlockParserLinearDoping import *
from .charonBlockParserAlbrechtMobility import *
from .charonBlockParserAugerRecombinationParameters import *
from .charonBlockParserDefectClusterParameters import *
from .charonBlockParserIntrinsicConc import *
from .charonBlockParserDoping import *
from .charonBlockParserUniformMoleFraction import *
from .charonBlockParserGaussDoping import *
from .charonBlockParserFarahmandMobility import *
from .charonBlockParserPhilipsMobility import *
from .charonBlockParserThermalConductivity import *
from .charonBlockParserPhilipsThomasMobility import *
from .charonBlockParserErfcDoping import *
from .charonBlockParserParticleStrikeParameters import *
from .charonBlockParserAvalancheGeneration import *
from .charonBlockParserGaAsMobility import *
from .charonBlockParserMOSFETMobility import *
from .charonBlockParserHeatCapacity import *
from .AroraMobility.AroraMobilityParserLib import *
from .GlobalMMSParameters.GlobalMMSParametersParserLib import *
from .StepDoping.StepDopingParserLib import *
from .MasettiMobility.MasettiMobilityParserLib import *
from .LucentMobility.LucentMobilityParserLib import *
from .HeatGeneration.HeatGenerationParserLib import *
from .CarrierLifetime.CarrierLifetimeParserLib import *
from .MGaussDoping.MGaussDopingParserLib import *
from .BandGapNarrowing.BandGapNarrowingParserLib import *
from .IncompleteIonizationDonor.IncompleteIonizationDonorParserLib import *
from .TrapSRH.TrapSRHParserLib import *
from .IncompleteIonizationAcceptor.IncompleteIonizationAcceptorParserLib import *
from .BulkFixedChargeParameters.BulkFixedChargeParametersParserLib import *
from .OpticalGeneration.OpticalGenerationParserLib import *
from .EffectiveDOS.EffectiveDOSParserLib import *
from .UniformDoping.UniformDopingParserLib import *
from .LinearDoping.LinearDopingParserLib import *
from .AlbrechtMobility.AlbrechtMobilityParserLib import *
from .AugerRecombinationParameters.AugerRecombinationParametersParserLib import *
from .DefectClusterParameters.DefectClusterParametersParserLib import *
from .IntrinsicConc.IntrinsicConcParserLib import *
from .Doping.DopingParserLib import *
from .UniformMoleFraction.UniformMoleFractionParserLib import *
from .GaussDoping.GaussDopingParserLib import *
from .FarahmandMobility.FarahmandMobilityParserLib import *
from .PhilipsMobility.PhilipsMobilityParserLib import *
from .ThermalConductivity.ThermalConductivityParserLib import *
from .PhilipsThomasMobility.PhilipsThomasMobilityParserLib import *
from .ErfcDoping.ErfcDopingParserLib import *
from .ParticleStrikeParameters.ParticleStrikeParametersParserLib import *
from .AvalancheGeneration.AvalancheGenerationParserLib import *
from .GaAsMobility.GaAsMobilityParserLib import *
from .MOSFETMobility.MOSFETMobilityParserLib import *
from .HeatCapacity.HeatCapacityParserLib import *



class MaterialBlockParserLib:
    "This is the  MaterialBlockParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "MaterialBlockParserLib"
        # create the linparser objects 
        self.lineParsers = []
        self.lineParsers.append(charonLineParserFixedMobility())
        self.lineParsers.append(charonLineParserAffinity())
        self.lineParsers.append(charonLineParserMaterialName())
        self.lineParsers.append(charonLineParserRadiativeRecombinationCoefficient())
        self.lineParsers.append(charonLineParserDopingFile())
        self.lineParsers.append(charonLineParserReferenceMaterial())
        self.lineParsers.append(charonLineParserIntrinsicConcentrationModel())
        self.lineParsers.append(charonLineParserRelativePermittivity())
        self.lineParsers.append(charonLineParserBandgap())
        self.lineParsers.append(charonLineParserMMSDoping())
        # create the blockparser objects 
        self.blockParsers = []
        self.blockParsers.append([charonBlockParserAroraMobility(),AroraMobilityParserLib()])
        self.blockParsers.append([charonBlockParserGlobalMMSParameters(),GlobalMMSParametersParserLib()])
        self.blockParsers.append([charonBlockParserStepDoping(),StepDopingParserLib()])
        self.blockParsers.append([charonBlockParserMasettiMobility(),MasettiMobilityParserLib()])
        self.blockParsers.append([charonBlockParserLucentMobility(),LucentMobilityParserLib()])
        self.blockParsers.append([charonBlockParserHeatGeneration(),HeatGenerationParserLib()])
        self.blockParsers.append([charonBlockParserCarrierLifetime(),CarrierLifetimeParserLib()])
        self.blockParsers.append([charonBlockParserMGaussDoping(),MGaussDopingParserLib()])
        self.blockParsers.append([charonBlockParserBandGapNarrowing(),BandGapNarrowingParserLib()])
        self.blockParsers.append([charonBlockParserIncompleteIonizationDonor(),IncompleteIonizationDonorParserLib()])
        self.blockParsers.append([charonBlockParserTrapSRH(),TrapSRHParserLib()])
        self.blockParsers.append([charonBlockParserIncompleteIonizationAcceptor(),IncompleteIonizationAcceptorParserLib()])
        self.blockParsers.append([charonBlockParserBulkFixedChargeParameters(),BulkFixedChargeParametersParserLib()])
        self.blockParsers.append([charonBlockParserOpticalGeneration(),OpticalGenerationParserLib()])
        self.blockParsers.append([charonBlockParserEffectiveDOS(),EffectiveDOSParserLib()])
        self.blockParsers.append([charonBlockParserUniformDoping(),UniformDopingParserLib()])
        self.blockParsers.append([charonBlockParserLinearDoping(),LinearDopingParserLib()])
        self.blockParsers.append([charonBlockParserAlbrechtMobility(),AlbrechtMobilityParserLib()])
        self.blockParsers.append([charonBlockParserAugerRecombinationParameters(),AugerRecombinationParametersParserLib()])
        self.blockParsers.append([charonBlockParserDefectClusterParameters(),DefectClusterParametersParserLib()])
        self.blockParsers.append([charonBlockParserIntrinsicConc(),IntrinsicConcParserLib()])
        self.blockParsers.append([charonBlockParserDoping(),DopingParserLib()])
        self.blockParsers.append([charonBlockParserUniformMoleFraction(),UniformMoleFractionParserLib()])
        self.blockParsers.append([charonBlockParserGaussDoping(),GaussDopingParserLib()])
        self.blockParsers.append([charonBlockParserFarahmandMobility(),FarahmandMobilityParserLib()])
        self.blockParsers.append([charonBlockParserPhilipsMobility(),PhilipsMobilityParserLib()])
        self.blockParsers.append([charonBlockParserThermalConductivity(),ThermalConductivityParserLib()])
        self.blockParsers.append([charonBlockParserPhilipsThomasMobility(),PhilipsThomasMobilityParserLib()])
        self.blockParsers.append([charonBlockParserErfcDoping(),ErfcDopingParserLib()])
        self.blockParsers.append([charonBlockParserParticleStrikeParameters(),ParticleStrikeParametersParserLib()])
        self.blockParsers.append([charonBlockParserAvalancheGeneration(),AvalancheGenerationParserLib()])
        self.blockParsers.append([charonBlockParserGaAsMobility(),GaAsMobilityParserLib()])
        self.blockParsers.append([charonBlockParserMOSFETMobility(),MOSFETMobilityParserLib()])
        self.blockParsers.append([charonBlockParserHeatCapacity(),HeatCapacityParserLib()])
        # create the parserLibrary objects 
        parserLibraries = []
        parserLibraries.append(AroraMobilityParserLib())
        parserLibraries.append(GlobalMMSParametersParserLib())
        parserLibraries.append(StepDopingParserLib())
        parserLibraries.append(MasettiMobilityParserLib())
        parserLibraries.append(LucentMobilityParserLib())
        parserLibraries.append(HeatGenerationParserLib())
        parserLibraries.append(CarrierLifetimeParserLib())
        parserLibraries.append(MGaussDopingParserLib())
        parserLibraries.append(BandGapNarrowingParserLib())
        parserLibraries.append(IncompleteIonizationDonorParserLib())
        parserLibraries.append(TrapSRHParserLib())
        parserLibraries.append(IncompleteIonizationAcceptorParserLib())
        parserLibraries.append(BulkFixedChargeParametersParserLib())
        parserLibraries.append(OpticalGenerationParserLib())
        parserLibraries.append(EffectiveDOSParserLib())
        parserLibraries.append(UniformDopingParserLib())
        parserLibraries.append(LinearDopingParserLib())
        parserLibraries.append(AlbrechtMobilityParserLib())
        parserLibraries.append(AugerRecombinationParametersParserLib())
        parserLibraries.append(DefectClusterParametersParserLib())
        parserLibraries.append(IntrinsicConcParserLib())
        parserLibraries.append(DopingParserLib())
        parserLibraries.append(UniformMoleFractionParserLib())
        parserLibraries.append(GaussDopingParserLib())
        parserLibraries.append(FarahmandMobilityParserLib())
        parserLibraries.append(PhilipsMobilityParserLib())
        parserLibraries.append(ThermalConductivityParserLib())
        parserLibraries.append(PhilipsThomasMobilityParserLib())
        parserLibraries.append(ErfcDopingParserLib())
        parserLibraries.append(ParticleStrikeParametersParserLib())
        parserLibraries.append(AvalancheGenerationParserLib())
        parserLibraries.append(GaAsMobilityParserLib())
        parserLibraries.append(MOSFETMobilityParserLib())
        parserLibraries.append(HeatCapacityParserLib())


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
