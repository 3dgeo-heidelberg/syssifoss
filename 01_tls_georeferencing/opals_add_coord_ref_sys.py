#!/bin/python
# -*- coding: utf-8 -*-
# Hannah Weiser (Heidelberg University)

"""
Script to add a coordinate reference system definition to the las file header.

example usage: opals_add_coord_ref_sys.py E:\SYSSIFOSS\pointclouds\*.laz 25832 E:\SYSSIFOSS\pointclouds_ref E:\SYSSIFOSS\04_opals_formatdef\LAS_1.4_TLS_ULS_amp_ref_dev.xml
"""

import os
import sys
import opals
from opals import Import, Export
import glob

def main(args):
    try:
        targ_file = args[0]
        epsg_code = args[1]
        result_dir = args[2]
        formatdef = args[3]

    except Exception as e:
        print("Usage: opals_add_coord_ref_sys.py target_file.las epsg_code result_dir\n"
              "target_file      is a las/laz file\n"
              "epsg_code        is the epsg code to add to the las header\n"
              
              "result_dir       the directory to which all interim and final results will be written\n"
              "formatdef        the path to an XML format definition file for Opals\n"
              "Execute script in the Opals Shell!\n"
              'example usage:"\n')
        print("Error: %s" % e)
        exit()

    infiles = []
    infiles.extend(glob.glob(targ_file))

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    epsg = "EPSG:%i" % int(epsg_code)

    for target_file in infiles:

        target_base = os.path.basename(target_file)
        target_fn = os.path.splitext(target_base)[0]

        imp = Import.Import()
        imp.inFile = target_file
        odm_file = os.path.join(result_dir, target_fn + ".odm")
        print(odm_file)
        imp.outFile = odm_file
        imp.tileSize = 0.5
        imp.iFormat = formatdef
        imp.run(reset=True)

        exp = Export.Export()
        exp.inFile = odm_file
        exp.outFile = os.path.join(result_dir, target_fn + ".laz")
        exp.globals.coord_ref_sys = epsg
        exp.oFormat = formatdef
        exp.tileSize = 0.5
        exp.run(reset=True)

        os.remove(odm_file)

if __name__ == '__main__':
    main(sys.argv[1:])