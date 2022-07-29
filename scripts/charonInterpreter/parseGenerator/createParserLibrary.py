

class createParserLibrary:

    def __init__(self):
        # This is a stateless class
       # First, create the parser class file
        self.nextLine = "\n"
        self.indent = "    "
        self.indent2 = self.indent+self.indent
        self.indent3 = self.indent2+self.indent
        self.indent4 = self.indent3+self.indent
        self.indent5 = self.indent4+self.indent
        self.indent6 = self.indent5+self.indent
        self.indent7 = self.indent6+self.indent
        self.indent8 = self.indent7+self.indent
        self.indent9 = self.indent8+self.indent
        self.indent10 = self.indent9+self.indent

    def writeParserLibraryClass(self,targetDirectory,libName,lineParsers,blockParsers):
        # create the init file so that the directory is treated as a module
        initFile = open(targetDirectory+"/__init__.py","w+")
        initFile.write(" ")
        initFile.close()

        #write the parserlibrary
        #First the name of the parserLibrary
        self.parserLibName = libName+"ParserLib"
        self.parserLibPath = targetDirectory+"/"
        parseLibFile = open(self.parserLibPath+self.parserLibName+".py","w+")

        fileContents = self.nextLine
        fileContents += "try:"+self.nextLine
        fileContents += self.indent+"import coloramaDISABLED as colors"+self.nextLine
        fileContents += "except ImportError:"+self.nextLine
        fileContents += self.indent+"class stubColors:"+self.nextLine
        fileContents += self.indent+"    \"subs for colors when colors doesn't exist on system\""+self.nextLine
        fileContents += self.indent+""+self.nextLine
        fileContents += self.indent+"    def __init__(self):"+self.nextLine
        fileContents += self.indent+"        self.Fore = colorClass()"+self.nextLine
        fileContents += self.indent+"        self.Back = colorClass()"+self.nextLine
        fileContents += self.indent+"        self.Style = styleClass()"+self.nextLine
        fileContents += self.indent+""+self.nextLine
        fileContents += self.indent+"class colorClass():"+self.nextLine
        fileContents += self.indent+"    \"stubbed color class\""+self.nextLine
        fileContents += self.indent+""+self.nextLine
        fileContents += self.indent+"    def __init__(self):"+self.nextLine
        fileContents += self.indent+"        self.BLACK = \"\""+self.nextLine
        fileContents += self.indent+"        self.BLUE = \"\""+self.nextLine
        fileContents += self.indent+"        self.WHITE = \"\""+self.nextLine
        fileContents += self.indent+"        self.RED = \"\""+self.nextLine
        fileContents += self.indent+"        self.GREEN = \"\""+self.nextLine
        fileContents += self.indent+""+self.nextLine
        fileContents += self.indent+"class styleClass():"+self.nextLine
        fileContents += self.indent+"    \"stubbed style class\""+self.nextLine
        fileContents += self.indent+""+self.nextLine
        fileContents += self.indent+"    def __init__(self):"+self.nextLine
        fileContents += self.indent+"        self.RESET_ALL = \"\""+self.nextLine
        fileContents += self.indent+"colors = stubColors()"+self.nextLine+self.nextLine
 

        fileContents += "import sys"+self.nextLine
        # import the line and block parser classes and any subordinate parser libraries
        for lP in lineParsers:
            fileContents += "from .charonLineParser"+lP+" import *"+self.nextLine
        # import the block parser classes
        for bP in blockParsers:
            fileContents += "from .charonBlockParser"+bP+" import *"+self.nextLine
        #import subordinate parser libraries
        for bP in blockParsers:
            fileContents += "from ."+bP+"."+bP+"ParserLib import *"+self.nextLine


        # Create the class
        fileContents += self.nextLine+self.nextLine+self.nextLine
        fileContents += "class "+self.parserLibName+":"+self.nextLine
        fileContents += self.indent+"\"This is the  "+self.parserLibName+" parser library \""+self.nextLine

        #Create Constructor
        fileContents += self.nextLine+self.nextLine
        fileContents += self.indent+"def __init__(self):"+self.nextLine
        #set the name
        fileContents += self.indent2+"# set the parser library name "+self.nextLine
        fileContents += self.indent2+"self.parserLibName = \""+self.parserLibName+"\""+self.nextLine
        #Create list of line parser objects
        fileContents += self.indent2+"# create the linparser objects "+self.nextLine
        fileContents += self.indent2+"self.lineParsers = []"+self.nextLine
        for pN in lineParsers:
            fileContents += self.indent2+"self.lineParsers.append(charonLineParser"+pN+"())"+self.nextLine

        #Create list of block parser and parser library objects
        fileContents += self.indent2+"# create the blockparser objects "+self.nextLine
        fileContents += self.indent2+"self.blockParsers = []"+self.nextLine
        for bN in blockParsers:
            fileContents += self.indent2+"self.blockParsers.append([charonBlockParser"+bN+"(),"+bN+"ParserLib()])"+self.nextLine

        #Create list of block parser objects
        fileContents += self.indent2+"# create the parserLibrary objects "+self.nextLine
        fileContents += self.indent2+"parserLibraries = []"+self.nextLine
        for bN in blockParsers:
            fileContents += self.indent2+"parserLibraries.append("+bN+"ParserLib())"+self.nextLine


        #Create a method that will loop through the line parsers and look for keyword match
        fileContents += self.nextLine+self.nextLine
        fileContents += self.indent+"def isThisMyLine(self,tokenizer,line):"+self.nextLine
        fileContents += self.indent2+"suggestedSyntax = []"+self.nextLine
        fileContents += self.indent2+"for lP in self.lineParsers:"+self.nextLine
        fileContents += self.indent3+"(self.isThisMe,suggestedSyntaxItem) = lP.isThisMe(tokenizer,line)"+self.nextLine
        fileContents += self.indent3+"if suggestedSyntaxItem is not \"\":"+self.nextLine
        fileContents += self.indent4+"suggestedSyntax.append(suggestedSyntaxItem)"+self.nextLine
        fileContents += self.indent3+"if self.isThisMe == True:"+self.nextLine
        fileContents += self.indent4+"return (True,lP,suggestedSyntax)"+self.nextLine
        fileContents += self.indent2+"return (False,None,suggestedSyntax)"+self.nextLine


        #Create a method that will loop through the block parsers and look for keyword match
        fileContents += self.nextLine+self.nextLine
        fileContents += self.indent+"def isThisMyBlock(self,tokenizer,line):"+self.nextLine
        fileContents += self.indent2+"for bP in self.blockParsers:"+self.nextLine
        fileContents += self.indent3+"self.isThisMe = bP[0].isThisMe(tokenizer,line)"+self.nextLine
        fileContents += self.indent3+"if self.isThisMe == True:"+self.nextLine
        fileContents += self.indent4+"return (True,bP[0],bP[1])"+self.nextLine
        fileContents += self.indent2+"return (False,None,None)"+self.nextLine


        #create a method that will loop through the library and generate help
        fileContents += self.nextLine+self.nextLine
        fileContents += self.indent+"def generateHelp(self,genHelp,indent):"+self.nextLine
        fileContents += self.indent2+"self.addIndent = \"     \""+self.nextLine
        fileContents += self.indent2+"cRStyle = \"\""+self.nextLine
        fileContents += self.indent2+"for lP in self.lineParsers:"+self.nextLine
        fileContents += self.indent3+"(self.helpLine,self.helpContent) = lP.getHelp(genHelp)"+self.nextLine
        fileContents += self.indent3+"self.helpContentList = self.helpContent.split(\"<>\")"+self.nextLine
        fileContents += self.indent3+"print (cRStyle+indent+colors.Fore.RED+colors.Back.WHITE+self.helpLine)"+self.nextLine
        fileContents += self.indent3+"cRStyle = \"\\n\""+self.nextLine
        fileContents += self.indent3+"for hCL in self.helpContentList:"+self.nextLine
        fileContents += self.indent4+"print (\"\\t\"+indent+colors.Fore.BLUE+colors.Back.WHITE+hCL.lstrip())"+self.nextLine
        fileContents += self.indent2+"for bP in range(len(self.blockParsers)):"+self.nextLine
        fileContents += self.indent3+"print (indent+colors.Fore.GREEN+colors.Back.WHITE+self.blockParsers[bP][0].getHelpLine().lstrip())"+self.nextLine
        fileContents += self.indent3+"self.blockParsers[bP][1].generateHelp(genHelp,indent+self.addIndent)"+self.nextLine
        fileContents += self.indent3+"print (indent+colors.Fore.GREEN+colors.Back.WHITE+self.blockParsers[bP][0].getHelpLine().replace(\"start\",\"end\").lstrip())"+self.nextLine

        fileContents += self.indent3+"print (indent+colors.Style.RESET_ALL)"+self.nextLine

        fileContents += self.nextLine+self.nextLine
 
        fileContents += self.indent+"def getName(self):"+self.nextLine
        fileContents += self.indent2+"return self.parserLibName"+self.nextLine

        parseLibFile.write(fileContents)

