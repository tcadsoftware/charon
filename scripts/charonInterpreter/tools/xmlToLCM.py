#! /usr/bin/env python3
from __future__ import print_function

import sys
import os
from os import path
from os import listdir

from .xmlToLCMConverter.xmlLCMConverter import *


xmlFilename = sys.argv[1]

x2lcm = xmlLCMConverter(xmlFilename)
x2lcm.convertFile()
x2lcm.writeLCMParameters()
