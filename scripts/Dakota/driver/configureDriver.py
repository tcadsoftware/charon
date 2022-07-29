
import copy
import sys

class configureDriver:
    "Configures the responses for the Charon driver for Dakota"


    def __init__(self):
        self.parametersDict = {}
        self.parametersDict['dakotaInputFilename'] = "NoName"
        self.parametersDict['executeProcs'] = 1
        self.responses = []
        self.selectedResponseObjects = []
        self.respArgList = []
        self.executeMethods = []
        self.calibrationConfig = []

    def configure(self):
        # Open the conf file
        configureInputList = list(open("driver.config"))
        # Parse the configuration to get filenames and responses
        for cIL in configureInputList:
            cILTokens = cIL.split()
            if len(cILTokens) == 0:  # Blank Line
                continue
            if cILTokens[0][0] == "#":  # Comment line
                continue

            if cILTokens[0].lower() == "executemethods":
                self.executeMethods.append(cILTokens[1:])

            if cILTokens[0].lower() == "dakotainput":
                self.parametersDict['dakotaInputFilename'] = cILTokens[1]

            if cILTokens[0].lower() == "executeprocs":
                self.parametersDict['executeProcs'] = cILTokens[1]

            if cILTokens[0].lower() == "response":
                self.responses.append(cILTokens[1])
                if len(cILTokens) > 2:  # arguments forthcoming
                    argList = cILTokens[2:]
                else:
                    argList = []
                self.respArgList.append(argList)

            if cILTokens[0].lower() == "calibration":
                self.calibrationConfig.append(cILTokens[1:])

        return (self.parametersDict,self.executeMethods,self.responses,self.respArgList,self.calibrationConfig)

    def setResponses(self,availableResponses):
        # Cycle through the available and responses list to return selected response object list

        for rL in self.responses:
            foundAvailable = False
            for aR in availableResponses:
                if aR.myResponse() == rL:
                    foundAvailable = True
                    self.selectedResponseObjects.append(copy.deepcopy(aR))
            if foundAvailable == False:
                print ("Error!! Cannot find an available response matching ",rL)
                print("Available Responses are:")
                for aR in availableResponses:
                    print (aR.myResponse())
                print ("Cannot Continue.  Exiting.")
                sys.exit(1)

        return self.selectedResponseObjects


    def validateAndOrderResponses(self,responses,inputParameters):

         ###################################################
        #  Read the expected responses from the input parameters passed from Dakota
        ###################################################
        inputs = list(open(inputParameters))
        dakResponses = []
        for iP in inputs:
            if "ASV_" in iP:
                dakResponses.append(iP.partition(":")[-1].strip())

        ###################################################
        #  Trim the Dakota responses for curve fits
        ###################################################

        curveTypes = []
        curveTypes.append("IVCurve")

        curvesPresent = []
        for index,rS in enumerate(responses):
            for cT in curveTypes:
                if rS.getMyType() == cT:
                    curvesPresent.append(rS.myResponse())
                    break

        modifiedDakResponses = []
        addedCurve = [False]*len(curvesPresent)
        for dK in dakResponses:
            foundLocal = False
            for index,cP in enumerate(curvesPresent):
                if dK.partition("_")[0] == cP and addedCurve[index] == False:
                    modifiedDakResponses.append(dK.partition("_")[0])
                    addedCurve[index] = True
                    foundLocal = True
                    break
            if foundLocal == False and dK.partition("_")[0] not in '  '.join(curvesPresent):
                modifiedDakResponses .append(dK)

        dakResponses = modifiedDakResponses

        foundResponse = [False]*len(dakResponses)
        returnResponses = [-1] * len(dakResponses)
        for index,dR in enumerate(dakResponses):
            for rs in responses:
                if dR.lower() == rs.myResponse().lower():
                    foundResponse[index] = True
                    returnResponses[index]  = rs

        for index,fR in enumerate(foundResponse):
            if fR == False:
                print ("Error!  Did not find ",dakResponses[index]," in the list of Charon Driver responses ",index)
                sys.exit(1)

        return returnResponses



    def readDakotaInputAndOrder(self,responses):
        filename = self.parametersDict['dakotaInputFilename']
        if filename == "NoName":  #Ignore this part and return an empty list
            return []
        processingResponses = False
        readNumExpected = False
        readResponses = False
        dakResponses = []
        numExpected = 0

        for resp in responses:
            print ("reordering ",resp.myResponse())

        dakotaIn = list(open(filename))
        for line in dakotaIn:
            lineTokens = line.split()
            if len(lineTokens) == 0:
                continue  ##empty line

            if lineTokens[0][0:9].lower() == "responses":  #Started responses block
                processingResponses = True

            if processingResponses == True:
                if lineTokens[0].lower() == "response_functions" or lineTokens[0].lower() == "calibration_terms":
                    numExpected = int(lineTokens[2])
                    readNumExpected = True
                if lineTokens[0].lower() == "descriptors" or lineTokens[0].lower() == "response_descriptors":
                    dakResponses.extend(lineTokens[2:])
                    readResponses = True

                if readResponses == True and readNumExpected == True:
                    processingResponses = False

        #Clean up responses
        for index,resp in enumerate(dakResponses):
            dakResponses[index] = resp[1:-1]

        #First order sanity check
        if len(dakResponses) != len(responses):
            print ("Error:  The number of responses specified in the driver is not equal to the number specified in the dakota input file.")
            print ("Error:  There are ",len(dakResponses)," in Dakota and ",len(responses)," specified in the driver.config")
            print ("Dakota Responses:")
            for resp in dakResponses:
                print (resp.myResponse())
            print ("Charon Driver responses:")
            for resp in responses:
                print (resp.myResponse())
            sys.exit(1)

        #Reorder selected response objects so that they correspond to the dakota ordering
        reorderIndex = [-1] * len(responses)
        for dRIndex,dakResp in enumerate(dakResponses):
            for RIndex,resp in enumerate(responses):
                if dakResp.lower() == resp.myResponse().lower():
                    reorderIndex[RIndex] = dRIndex

        for index,rOI in enumerate(reorderIndex):
            if rOI < 0:
                print ("Error!  Did not find ",responses[index].myResponse()," in the list of Dakota responses ",index)
                sys.exit(1)

        returnResponses = [-1] * len(responses)
        for index,rI in enumerate(reorderIndex):
            returnResponses[rI] = responses[index]

        return returnResponses



