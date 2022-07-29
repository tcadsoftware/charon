#! /usr/bin/env python3



# This script compares dakota output for use in testing.


import os
import sys

verbose = False

if len(sys.argv) != 3:
    print("Error!  This script MUST take two arguments.  compareDakotaOutput.py file1 file2")
    sys.exit(1)

filename1 = sys.argv[1]
filename2 = sys.argv[2]

file1Lines = list(open(filename1))
file2Lines = list(open(filename2))

if len(file1Lines) != len(file2Lines):
    print ("Error! the files are different sizes.  Comparison is not possible")
    sys.exit(1)

if len(file1Lines[0].split()) != len(file2Lines[0].split()):
    print ("Error! the files have a different number of columns.  Comparison is not possible")
    sys.exit(1)


#The first line of the file is column headers.  Check to be sure they're the same
errorList = []
for index,f1 in enumerate(file1Lines[0].split()):
    if f1 != file2Lines[0].split()[index]:
        errorList.append([index,f1,file2Lines[0].split()[index]])


if len(errorList) != 0:
    print("There are discrepancies in the dakota output column headers.  Most likely this means that the Dakota run has changed and there is no hope that the output will match.")
    print("---------------------------------------------------------------------------------")
    for err in errorList:
        print ("In column ",err[0],err[1]," != ",err[2])


defaultFloor = 1e-12
defaultRelTol = 1e-10


columnNameList = file2Lines[0].split()[2:]
floor = []
relTol = []

for col in columnNameList:
    if verbose == True:
        print ("Comparing ",col)
    floor.append(defaultFloor)
    relTol.append(defaultRelTol)

######################################

errorList.clear()
offset =  2

for lineIndex,line in enumerate(file1Lines[1:]):
    f1 = line.split()
    f2 = file2Lines[lineIndex+1].split()
    if verbose == True:
        print ("Comparing LInes:")
        print ("     ",f1)
        print ("     ",f2)
    if verbose == True:
        print ("Comparing lines ",lineIndex)
    for index,f1vs in enumerate(f1[offset:]):
        f1v = float(f1vs)
        f2v = float(f2[index+offset])
        if verbose == True:
            print ("Comparing ",f1v,f2v)
        if f1v < floor[index] and f2v < floor[index]:
            continue
        relError = (f1v-f2v)/f1v
        if  relError > relTol[index]:
            errorList.append( [lineIndex+1, index+offset, columnNameList[index], f1v, f2v, relError])

######################################


if len(errorList) != 0:
    print ("There are errors in the Dakota output.")
    for err in errorList:
        print ("At [",err[0],",",err[1],"] ",err[2]," ~ ",err[3]," , ",err[4]," - relative error = ",err[5],sep='')

    print ("Dakota comparison failed")

    sys.exit(1)

else:

    print("files are the same")

