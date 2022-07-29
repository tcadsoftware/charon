
import sys

from modules.commonFunctions.parse.parse import *

class parseExecuteMethod:
    "parse the execute method line for preprocessing and execution"

    def __init__(self):
        self.executeMethod = []


    def parseEM(self,executeMethod):
        self.executeMethod = executeMethod

        
        app = ""
        template = ""
        options = ""
        # Do some processing on the list

        args = " ".join(self.executeMethod)
        args = " ".join(args.replace("="," = ").split())        
        args = " ".join(args.replace(", ",",").split())        
        args = " ".join(args.replace(" ,",",").split())        

        appD = search("app = {app:>S}",args)
        templateD = search("template = {template:>S}",args)
        optionsD = search("options = {options:>S}",args)

        if appD == None:
            print ("Error! No application specified in the executeMethod line.")
            sys.exit(1)
        else:
            app = appD['app']

        if templateD == None:
            print ("Error! No template file specified in the executeMethod line.")
            sys.exit(1)
        else:
            template = templateD['template']
            if "," in template:
                template = template.replace(",","  ")

            template = template.split()

        if optionsD != None:
            options = "--"+(optionsD['options']).replace(","," --")


        return (app,template,options)
