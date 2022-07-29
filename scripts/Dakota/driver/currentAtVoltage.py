#! /usr/bin/env python3


import os
import sys
import argparse

from modules.currentAtVoltage.currentAtVoltage import *


#######################################################
# Create the argument parser
#######################################################

parser = argparse.ArgumentParser()

#######################################################
##  add command line options 
#######################################################

parser.add_argument("-v","--voltage",help="Voltage at which to compute current")
parser.add_argument("-f","--filename",help="File name that holds I-V data to compute current.")
parser.add_argument("-vc","--voltColumn",help="Specify the column for the voltage sweep starting from index 0")
parser.add_argument("-cc","--currentColumn",help="Specify the column for the current starting from index 0")
args = parser.parse_args()


#######################################################
##  process command line options 
#######################################################

targetVoltage = -1
filenameSpecified = False
filename = "NoName"
voltageSpecified = False
voltColumnSpecified = False
currentColumnSpecified = False
voltColumn = ""
currentColumn = ""


if args.voltage:
    targetVoltage = args.voltage
    voltageSpecified = True

if args.filename:
    filenameSpecified = True
    filename = args.filename

if args.voltColumn != None:
    voltColumnSpecified = True
    voltColumn = args.voltColumn

if args.currentColumn != None:
    currentColumnSpecified = True
    currentColumn = args.currentColumn


#######################################################
##  sanity checks
#######################################################

if not voltageSpecified:
    print ("Error! You must specify a target voltage.")
    sys.exit(1)

if not filenameSpecified:
    print ("Error!  You must specify a file name for the I-V data!")
    sys.exit(1)


cAV = currentAtVoltage()

respArgs = []

voltageArg = "voltage = "+str(targetVoltage)
respArgs.append(voltageArg)

filenameArg = " filename="+filename
respArgs.append(filenameArg)

if voltColumnSpecified == True:
    voltColumnArg = " voltColumn = "+voltColumn
    respArgs.append(voltColumnArg)

if currentColumnSpecified == True:
    currentColumnArg = " currentColumn = "+currentColumn
    respArgs.append(currentColumnArg)

cAV.createDict(respArgs)

current = cAV.getResponse()

print ("current = ",current)


