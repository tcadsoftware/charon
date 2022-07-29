#! /usr/bin/env python3

import glob
import shutil

file = glob.glob('meshRefine*.so')
shutil.move(file[0],'../meshRefine.so')

