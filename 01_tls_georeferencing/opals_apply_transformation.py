#!/bin/python
# -*- coding: utf-8 -*-
# Hannah Weiser (Heidelberg University)

"""
Script to apply several transformations, specified in a comma-separated string to one las file

example usage: opals_apply_transformation.py "E:\SYSSIFOSS\MODEL_TREES\F_BR08_01.las" "Affine[-477800.00 -5428900.00],Affine[0.996834 0.079511 -0.079511 0.996834],Affine[477790.775785 5428916.817188 -10.56571],affine[0.999993326765 0.003648593929 0.000184900030 -0.003648699619 0.999993178431 0.000574527986 -0.000182802549 -0.000575198797 0.999999817865 -19804.927977241459 1780.504732256755 3210.432627653161]" "E:\SYSSIFOSS\06_opals_formatdef\LAS_1.4_TLS_ULS_amp_ref_dev.xml"
Important: adapt the variable formatdef to point to the format definition .xml on your local computer

Usage (in Opals Shell):
python opals_apply_transformation.py <target_file.las> <trafo_string> <result_dir> <format_def.xml>
    target_file      is a las file (from TLS point cloud)
    trafo_string     is opals filter strings for multiple translations and rotations, each separated by comma,
                        e.g., "Affine[-477800.00 -5428900.00],Affine[0.996834 0.079511 -0.079511 0.996834],Affine[477790.775785 5428916.817188 -10.56571],affine[0.999993326765 0.003648593929 0.000184900030 -0.003648699619 0.999993178431 0.000574527986 -0.000182802549 -0.000575198797 0.999999817865 -19804.927977241459 1780.504732256755 3210.432627653161]"
    result_dir       is the directory to which all interim and final results will be written
    formatdef        is an XML format definition file for Opals,
                         e.g., E:\SYSSIFOSS\06_opals_formatdef\LAS_1.4_TLS_ULS_amp_ref_dev.xml
"""

import os, sys
import opals
from opals import Export, Import
import glob

def main(args):
    try:
        targ_file = args[0]
        trafo_string = args[1]
        result_dir = args[2]
        formatdef = args[3]

    except Exception as e:
        print("Usage: <target_file.las> <trafo_string> <result_dir> <format_def.xml>\n"
              "target_file      the path to the las file (from TLS point cloud)\n"
              "trafo_string     opals filter strings for translations and rotations separated by comma\n"
              "result_dir       the directory to which all interim and final results will be written\n"
              "formatdef        the path to an XML format definition file for Opals\n"
              "Execute script in the Opals Shell!\n"
              'example usage: ApplyTransformation.py "E:\SYSSIFOSS\MODEL_TREES\F_BR08_01.las" "Affine[-477800.00 -5428900.00],Affine[0.996834 0.079511 -0.079511 0.996834],Affine[477790.775785 5428916.817188 -10.56571],affine[0.999993326765 0.003648593929 0.000184900030 -0.003648699619 0.999993178431 0.000574527986 -0.000182802549 -0.000575198797 0.999999817865 -19804.927977241459 1780.504732256755 3210.432627653161]" "E:\SYSSIFOSS\REGISTRATION\F_BR08" "E:\SYSSIFOSS\06_opals_formatdef\LAS_1.4_TLS_ULS_amp_ref_dev.xml"\n')
        print("Error: %s" % e)
        exit()

    infiles = []
    infiles.extend(glob.glob(targ_file))

    trafo_string = trafo_string.split(',')
    print(trafo_string)

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    for target_file in infiles:

        target_base = os.path.basename(target_file)
        target_fn = os.path.splitext(target_base)[0]

        imp = Import.Import()
        imp.inFile = target_file
        target_current = os.path.join(result_dir, target_fn + "transformed1.odm")
        imp.outFile = target_current
        imp.filter = "%s" % trafo_string[0]
        imp.tileSize = 0.5
        imp.iFormat = formatdef
        imp.run(reset=True)

        to_remove = []
        to_remove.append(target_current)
        for i, trafo in enumerate(trafo_string):
            if i == 0:
                continue
            elif i == len(trafo_string)-1:
                exp = Export.Export()
                exp.inFile = target_current
                target_current = os.path.join(result_dir, target_fn + "_t.las")
                exp.outFile = target_current
                exp.filter = "%s" % trafo
                exp.oFormat = formatdef
                exp.tileSize = 0.5
                exp.run(reset=True)
            else:
                exp = Export.Export()
                exp.inFile = target_current
                target_current = os.path.join(result_dir, target_fn + "transformed%s.odm" % (i+1))
                to_remove.append(target_current)
                exp.outFile = target_current
                exp.filter = "%s" % trafo
                exp.globals.coord_ref_sys = "EPSG:25832"
                exp.tileSize = 0.5
                exp.run(reset=True)
            os.remove(to_remove[i-1])

if __name__ == '__main__':
    main(sys.argv[1:])