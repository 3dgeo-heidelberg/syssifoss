#!/bin/python
# -*- coding: utf-8 -*-
# Hannah Weiser (Heidelberg University)
# h.weiser@stud.uni-heidelberg.de


import os, sys
import glob
import opals
from opals import Export, Import, ICP

def main(args):
    try:
        fixFile = args[0]
        movFile_ds = args[1]
        result_dir = args[2]

    except Exception as e:
        print("Usage: opals_run_icp.py fixFile.odm movFile_ds.odm result_dir\n "
              "fixFile is an odm file (from ULS point cloud)\n"
              "movFile_ds is an odm file (downsampled and coarsely registered TLS point cloud)\n"
              "all interim and final results will be stored in result_dir\n"
              "execute script in Opals Shell in the directory containing the input files\n"
              "use wildcards for movFile_ds to select several mov_files belonging to the same fixFile\n"
              "example: opals_run_icp.py ..\ULS_BR03_merged.odm ..\*_P03_*.odm ...\coreg")
        print("Error: %s" % e)
        exit()

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    mov_files = []
    mov_files.extend(glob.glob(movFile_ds))

    for movFile in mov_files:
        mov_file_base = os.path.basename(movFile)
        mov_file_name = os.path.splitext(mov_file_base)[0]

        icp = ICP.ICP()
        icp.inFile = fixFile, movFile
        icp.outDirectory = result_dir + "\icpOut_" + mov_file_name
        icp.tempDirectory = result_dir + "\icpTemp_" + mov_file_name
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
