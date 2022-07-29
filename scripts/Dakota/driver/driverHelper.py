#! /usr/bin/env python3


class driverHelper:
    "Helper class for response functions"


    def __init__(self,responses):
        self.responses = responses


    def getHelp(self,request):
        foundHelp = False
        for resp in self.responses:
            if resp.myResponse().lower() == request.lower():
                resp.myResponse()
                resp.myHelp()
                foundHelp = True

        if foundHelp == False:
            print ("Error!  I did not find help for ",request)



