#! /usr/bin/env python3


import os
import sys
import argparse

from modules.thresholdVoltage.thresholdVoltage import *


#######################################################
# Create the argument parser
#######################################################

parser = argparse.ArgumentParser()

#######################################################
##  add command line options 
#######################################################

#parser.add_argument("filename",help="Specify the intpreter input file name.")
parser.add_argument("-m","--method",help="method to use to compute Vt (slope, current)")
parser.add_argument("-c","--current",help="current at which to compute Vt if method is current.")
parser.add_argument("-f","--filename",help="File name that holds I-V data to compute Vt.")
parser.add_argument("--icap",help="If method is slope, it will limit the evaluation to a current no higher than specified")
parser.add_argument("-vc","--voltColumn",help="Specify the column for the voltage sweep starting from index 0")
parser.add_argument("-cc","--currentColumn",help="Specify the column for the current starting from index 0")

args = parser.parse_args()


#######################################################
##  process command line options 
#######################################################

useCurrent = False
useSlope = False
targetCurrent = -1
filenameSpecified = False
filename = "NoName"
currentSpecified = False
method = ""
icapSpecified = False
voltColumnSpecified = False
currentColumnSpecified = False
voltColumn = ""
currentColumn = ""

if args.method != None:
    method = args.method.lower()
    if method.lower() == "slope":
        useSlope = True
    elif method.lower() == "current":
        useCurrent = True
else:
    print("Error!  You must specify a method!")
    sys.exit(1)


if args.current:
    targetCurrent = args.current
    currentSpecified = True

if args.filename:
    filenameSpecified = True
    filename = args.filename

if args.icap:
    icapSpecified = True
    icap = args.icap

if args.voltColumn != None:
    voltColumnSpecified = True
    voltColumn = args.voltColumn

if args.currentColumn != None:
    currentColumnSpecified = True
    currentColumn = args.currentColumn




#######################################################
##  sanity checks
#######################################################

if method.lower() == "current" and not useCurrent:
    print ("Error!  If you use the current method, you must specify a current!")
    sys.exit(1)


if not filenameSpecified:
    print ("Error!  You must specify a file name for the I-V data!")
    sys.exit(1)


Vt = thresholdVoltage()

respArgs = []

filenameArg = " filename="+filename

if method == "slope":
    methodArg = " method = slope"
    respArgs.append(methodArg)

if method == "current":
    currentArg = " current = "+targetCurrent
    respArgs.append(currentArg)

# Set standalone
standaloneArg = "standalone = True"
respArgs.append(standaloneArg)

filenameArg = " filename="+filename
respArgs.append(filenameArg)

if icapSpecified == True:
    icapSpecArg = "icap = "+icap
    respArgs.append(icapSpecArg)

if voltColumnSpecified == True:
    voltColumnArg = " voltColumn = "+voltColumn
    respArgs.append(voltColumnArg)

if currentColumnSpecified == True:
    currentColumnArg = " currentColumn = "+currentColumn
    respArgs.append(currentColumnArg)



Vt.createDict(respArgs)

thresh = Vt.getResponse()

threshTokens = thresh.split()

#print the results
print ("Vt = ",threshTokens[0]," at ",threshTokens[1]," amps")



