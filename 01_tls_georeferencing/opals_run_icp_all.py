#!/bin/python
# -*- coding: utf-8 -*-
# Hannah Weiser (Heidelberg University)
# h.weiser@stud.uni-heidelberg.de

"""
Script to run Opals ICP with all input files (given using wildcards)
Each movFiles is matched with the correct fixFile using the plot ID, contained in the filename (lines 49 to 52)
The index in line 49 has to be adapted, if the ID is expected at a different location within the filename
"""

import os, sys
import glob
import opals
from opals import Export, Import, ICP

def main(args):
    try:
        fix_file = args[0]
        mov_file_ds = args[1]
        result_dir = args[2]

    except Exception as e:
        print("Usage: opals_run_icp_all.py fix_file.odm mov_file_ds.odm result_dir\n "
              "fix_file is an odm file (from ULS point cloud)\n"
              "mov_file_ds is an odm file (downsampled and coarsely registered TLS point cloud)\n"
              "all interim and final results will be stored in result_dir\n"
              "execute script in Opals Shell in the directory containing the input files\n"
              "use wildcards to select multiple input files\n"
              "example: opals_run_icp_all.py dir\ULS_*.odm dir\2019*.odm dir")
        print("Error: %s" % e)
        exit()

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    fixFiles = []
    fixFiles.extend(glob.glob(fix_file))

    movFiles = []
    movFiles.extend(glob.glob(mov_file_ds))

    for movFile in movFiles:
        movFile_base = os.path.basename(movFile)
        movFile_name = os.path.splitext(movFile_base)[0]

        ID = movFile_name[10:12]

        matching = [path for path in fixFiles if ID in path]
        fix_file = matching[0]

        icp = ICP.ICP()
        icp.inFile = fix_file, movFile
        icp.outDirectory = result_dir + "\icpOut_" + movFile_name
        icp.tempDirectory = result_dir + "\icpTemp_" + movFile_name
        icp.trafoType = opals.Types.LSMTrafoType.rigid
        icp.voxelSize = 1
        icp.samplingDist = 0.5
        icp.searchRadius = 0.2
        icp.maxAngleDev = 5
        icp.maxIter = 5
        icp.maxSigma = 0.1
        icp.transformData = 0
        icp.run()
        icp.reset()

if __name__ == '__main__':
    main(sys.argv[1:])
