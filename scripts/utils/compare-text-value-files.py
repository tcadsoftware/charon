#!/usr/bin/env python
import sys
import re
import csv

# Script to compare two text files containing columnar floating point data.  In
# Charon this is currently used to compare two files containing Gummel data.
# The columns must be space-delimited, and the first row must contain the
# titles of the columns.  Optional column demarcation characters (e.g. |) are
# permitted.

DEBUG=0

try:
  if not (len(sys.argv) == 3 or len(sys.argv) == 4):
    raise Exception("Failure: Command line - " + sys.argv[0] + " [filename 1] [filename 2] <tolerance>")

  if (len(sys.argv) == 4):
    tolerance = float(sys.argv[3])
  else:
    tolerance = float(1.0e-6)

  if (DEBUG != 0):
    print("[DBG]:tolerance = " + str(tolerance))

  # Read the first file
  file1_row_count=0
  file1_values = []
  with open(sys.argv[1], 'r') as csvfile:
    filereader = csv.reader(csvfile, delimiter=' ')
    for row in filereader:
      if (file1_row_count == 0):
        file1_titles = row
        if (DEBUG != 0):
          countit=0
          for entry in file1_titles:
            print("[DBG]: file 1: Title column " + str(countit) + ": " + entry)
            countit = countit + 1
      else:
        file1_values.append(row)
      file1_row_count = file1_row_count + 1

  # Read the second file
  file2_row_count=0
  file2_values = []
  with open(sys.argv[2], 'r') as csvfile:
    filereader = csv.reader(csvfile, delimiter=' ')
    for row in filereader:
      if (file2_row_count == 0):
        file2_titles = row
        if (DEBUG != 0):
          countit=0
          for entry in file2_titles:
            print("[DBG]: file 2: Title column " + str(countit) + ": " + entry)
            countit = countit + 1
      else:
        file2_values.append(row)
      file2_row_count = file2_row_count + 1


  # If number of entries in the files doesn't match, files aren't a match
  if (file1_row_count != file2_row_count):
    raise Exception("Failure: files have a different number of rows")

  # Compare the column titles
  for i in range(len(file1_titles)):
    if (file1_titles[i] != file2_titles[i]):
      raise Exception("Column " + str(i) + " titles do not match. " + file1_titles[i] + ":" + file2_titles[i])

  # Compare the values. IF they're non-zero then a relative difference
  # is performed to the requested tolerance. Otherwise an absolute
  # diff is done.
  for i in range(len(file1_values)):
    if (len(file1_values[i]) != len(file2_values[i])):
      raise Exception("Row " + str(i) + " contains a different number of columns between the two files")
    
    for j in range(len(file1_values[i])):
      if (file1_values[i][j] == file2_values[i][j]):
        continue
      val1 = float(file1_values[i][j])
      val2 = float(file2_values[i][j])
      diffval = val1 - val2
      if (val1 != float(0.0)):
        diffval = diffval / val1

      diffval = abs(diffval)
      if (diffval > tolerance):
        raise Exception("Entries [" + str(file1_values[i][j]) + ", " +
                        str(file2_values[i][j]) + "] at row:" + str(i) + " col:" + str(j) + " differ")

except Exception as e:
  print("Test Failed: ")
  print(e)
else:
  print("Test Passed")
  
