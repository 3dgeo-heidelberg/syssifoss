#!/bin/python
# -*- coding: utf-8 -*-
# Hannah Weiser (Heidelberg University), July 2020
# h.weiser@stud.uni-heidelberg.de


"""
This scripts iterates through LAS files of tree point clouds and computes the tree position.

The position is derived as:

x,y coordinates: center of the bounding rectangle of the stem slice from the lowest point to 0.5 m
z coordinate: DTM height value of the pixel at the derived xy-position
Furthermore, coordinates are transformed from the specified crs (default: EPSG:25832) to WGS84

Usage: python tree_positions.py <data\*.laz> <grids\*.tif> <out_dir>
    inf             LAS/LAZ file(s) (supports wildcards)
    inf_raster      DTM(s) in geotiff format
                    - the same tree ID has to be contained in the filename of infile.las and translation.txt
                    - the tree ID is expected in the format: speciesID_plotID_number
    out_dir         directory to write positions to (one file per tree)

Example usage:
"""


import numpy as np
import os
import sys
from laspy.file import File
import glob
import rasterio
from pyproj import Transformer
import math
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist, squareform

inf = sys.argv[1]
inf_raster = sys.argv[2]
out_dir = sys.argv[3]
if len(sys.argv) > 4:
    epsg = sys.argv[4]
else:
    epsg = 25832

infiles = glob.glob(inf)
rasterfiles = glob.glob(inf_raster)

transformer = Transformer.from_crs(epsg, 4326)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

for i, fid in enumerate(infiles):
    if i % 10 == 0:
        if i == 0:
            print("Start processing...")
        else:
            print("Processing point cloud %d..." % i)
    file, ext = os.path.splitext(os.path.basename(fid))
    fileparts = file.split("_")
    plot_ID = fileparts[1]

    # get corresponding raster file (based on ID)
    raster_file = [r for r in rasterfiles if plot_ID in r]

    out = os.path.join(out_dir, file.replace("_20", "_m_20") + "_position_epsg%s.txt" % epsg)
    outfile = open(out, "w")
    inFile = File(fid, mode="r")
    point_cloud = np.array([inFile.x, inFile.y, inFile.z])
    pts = point_cloud.T
    minz = np.min(pts[:, 2])

    h = 1.0  # 2.0
    dZ = 0.2/2
    x = None
    while x is None:
        stem_slice = pts[(pts[:, 2] > (minz+(h-dZ))) & (pts[:, 2] < (minz+h+dZ))]
        if stem_slice.shape[0] > 2:
            x = np.sum(stem_slice[:, 0]/stem_slice.shape[0])
            y = np.sum(stem_slice[:, 1]/stem_slice.shape[0])

        h += dZ*2

    with rasterio.open(raster_file[0]) as src:
        z = list(src.sample([[x, y]]))
        z = z[0]

    outfile.write("%.6f,%.6f,%.6f\n" % (x, y, z))
    outfile.close()

    x, y, z = transformer.transform(x, y, z)

    if not os.path.exists(out_dir.replace("epsg%s" % epsg, "WGS84")):
        os.makedirs(out_dir.replace("epsg%s" % epsg, "WGS84"))

    out = os.path.join(out_dir.replace("epsg%s" % epsg, "WGS84"), file + "_position.txt")
    outfile = open(out, "w")
    outfile.write("%.6f,%.6f,%.6f\n" % (x, y, z))

    outfile.close()
    inFile.close()
