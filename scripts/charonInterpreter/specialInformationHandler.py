
import os
import sys

class specialInformationHandler:
    "A very simple class to extract special information from a list"


    def __init__(self):
        self.key = ""
        self.returnValue = []

    def getValue(self,specialInformation,keySearch):
        for keyList in specialInformation:
            if keyList[0] == keySearch:
                return keyList[1:]

            print ("ERROR!  No keys in special informatin match ",keySearch)


    def cleanTextData(self,specialInformation):
        self.locaDefault = "currents-loca.dat"
        self.transientDefault = "current_time.csv"
        deleteFilename = ""
        for keyList in specialInformation:
            if keyList[0] == "textDataFilename":
                deleteFilename = keyList[1]

        if deleteFilename != "":
            if os.path.isfile(deleteFilename):
                os.remove(deleteFilename)
        else: # delete by default name
            if os.path.isfile(self.locaDefault):
                os.remove(self.locaDefault)
            if os.path.isfile(self.transientDefault):
                os.remove(self.transientDefault)

