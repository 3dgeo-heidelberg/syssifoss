#!/bin/python
# -*- coding:utf8 -*-
# M. Bruggisser (TU Vienna), H. Weiser (Heidelberg University)
# h.weiser@stud.uni-heidelberg.de

"""
This program uses Opals (https://opals.geo.tuwien.ac.at/html/stable/ref_parameters.html).
The script extracts individual stem positions from a UAV-borne laser scanning (ULS) point cloud.
Steps:
1. Import of ULS LAS or LAZ file
2. Computation of normals and linearity: linearity = sqrt(1-NormalEigenvalue2/NormalEigenvalue1)
3. Filtering by linearity and z component of normal
4. Computation of point count and filtering by point count
5. Stem segmentation
6. Extraction of stem position (2D center of gravity of stem segment) and Export as .txt-file

Note: There are different ways to tweak the script to increase or decrease the sensitivity and adapt it to specific
input data or specific aims (i.e. extracting only well sampled stems etc.)
a) adapting the point count threshold (line 83)
b) adapting the minimum segment size (line 89)
c) adapting the normalZ threshold (esp. if also interested in tilted stems) and maybe the linearity threshold (line 68)

Execution (in OPALS Shell):
python opals_stem_detection_uls.py <input.las> <path to folder> <stemCenters.xyz>

<input.las> is the path to the input point cloud from which to extract the stem positions
<path to folder> is the working directory, in which intermediate and final results will be written
<stemCenters.xyz> is the name of the output txt-file, to which stem positions will be written
    -> one line per stem, x and y position, separated by comma
"""

import opals
from opals import Cell, Grid, AddInfo, Normals, Export, Algebra, PointStats, Import, PointStats, Segmentation, pyDM
import os, sys
import numpy as np

def main(base_file, work_dir, outID):
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    # import las file
    imp = Import.Import()
    imp.inFile = base_file
    base_file = base_file.replace('.las', '.odm').replace('laz', '.odm')
    imp.outFile = os.path.join(work_dir, base_file)
    imp.run(reset=True)

    # compute normals
    nor = Normals.Normals()
    nor.inFile = base_file
    nor.searchRadius = 1
    nor.neighbours = 9999
    nor.storeMetaInfo = opals.Types.NormalsMetaInfo.medium
    nor.searchMode = opals.Types.SearchMode.d3
    nor.run(reset=True)

    # compute linearity
    add = AddInfo.AddInfo()
    add.inFile = base_file
    add.globals.points_in_memory = 1000000
    add.attribute = "_linearity(float)=sqrt(1-NormalEigenvalue2/NormalEigenvalue1)"
    add.run(reset=True)

    # export stem points by filtering with the variables _linearity and NormalZ
    exp = Export.Export()
    exp.inFile = base_file
    exp.outFile = os.path.join(work_dir, "slice_stems.odm")
    exp.filter = "Generic[_linearity > 0.92] and Generic[NormalZ < 0.13]"
    exp.run(reset=True)

    # calculate point count
    ps = PointStats.PointStats()
    ps.inFile = os.path.join(work_dir, "slice_stems.odm")
    ps.searchRadius = 0.5
    ps.searchMode = opals.Types.SearchMode.d2
    ps.feature = "pcount"
    ps.run(reset=True)

    # export final stem points by applying a point count threshold
    exp = Export.Export()
    exp.inFile = os.path.join(work_dir, "slice_stems.odm")
    exp.outFile = os.path.join(work_dir, "slice_stems_over10.odm")
    exp.filter = "Generic[_dZPcount > 10]"
    exp.run(reset=True)

    # perform stem segmentation
    seg = Segmentation.Segmentation()
    seg.inFile = os.path.join(work_dir, "slice_stems_over10.odm")
    seg.minSegSize = 10
    seg.searchRadius = 0.5
    seg.run(reset=True)


    # open the odm / pyDM.Datamanager.load parameters: filename(string), readOnly(bool) threadSafety(bool)
    dm = pyDM.Datamanager.load(os.path.join(work_dir, "slice_stems_over10.odm"), True, False)

    psize = dm.sizePoint()

    # [x,y,z,id]
    pc = np.zeros([psize,4])

    # create an attribute layout for iteration (optional)
    lf = pyDM.AddInfoLayoutFactory()
    lf.addColumn(pyDM.ColumnSemantic.SegmentID)
    layout = lf.getLayout()

    # now iterate over all points
    for idx, pt in enumerate(dm.points(layout)):
        # output point
        pc[idx,0] = pt.x
        pc[idx,1] = pt.y
        pc[idx,2] = pt.z
        pc[idx,3] = pt.info().get(0)

    lbl = np.unique(pc[:,-1])
    lbl = lbl[np.logical_not(np.isnan(lbl))]

    ulsPC = np.zeros([len(lbl), 2], dtype = float)

    for ii,ll in enumerate(lbl):
        pts = pc[pc[:,3] == ll, 0:2]
        ulsPC[ii, 0] = np.mean(pts[:, 0]) #mean x
        ulsPC[ii, 1] = np.mean(pts[:, 1]) #mean y

    outID = os.path.join(work_dir,outID)
    np.savetxt(outID, ulsPC, delimiter=',')

if __name__=='__main__':
    main(*sys.argv[1:])  # <input.las> <path to folder> <stemCenters.xyz>