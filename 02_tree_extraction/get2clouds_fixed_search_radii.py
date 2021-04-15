#!/bin/python
# -*- coding: utf-8 -*-

"""
Grab XYZ point attributes from source point cloud using points from target point cloud:
  - For each point in the target point cloud, the closest k points in the source point cloud within a specified
    search_radius are written to the output point cloud (with attributes from source point cloud)
  - Assignment is done via Euclidean distance : fixed search radii are used: 0.1 m, 0.3 m, 0.5 m amd 1.0 m
        -> written to separate folders

Example use case: Extract single tree point cloud from an ALS plot point cloud using high resolution ULS or TLS
target point clouds of single trees

Important:
- One source point cloud and multiple target point clouds can be provided.
- Point clouds are expected in LAS or LAZ format
  (laszip.dll, e.g. from LAStools has to be added to the path to read LAZ files)
- If source point cloud is LAZ, output will also be LAZ

(c) H Weiser, 2020 from Bernhard Hoefle, 2019
"""

import sys
import argparse
from pykdtree.kdtree import KDTree      # Install: https://github.com/storpipfugl/pykdtree
import numpy as np
import pandas as pd                     # pandas needs to be installed: pip install pandas
import time
from laspy.file import File             # laspy needs to be installed: pip install laspy
import math
import os
import glob


def readLas(file):
    source_cloud = np.array([file.x, file.y, file.z]).T
    return source_cloud


def readUserArguments():
    # Generate help and user argument parsing
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="(c) Bernhard Hoefle (2019) - 3DGeo @ Heidelberg University")
    parser.add_argument("--source", dest='source', type=str, help="Source point cloud from which attributes and points shall be taken. It can be in ascii or las format", required=True)
    parser.add_argument("--target", dest='target', nargs='+', type=str, help="Target point cloud(s) to which the attributes from source point cloud should be appended and which is the spatial select: only points of source close to points of the target are considered", required=True)
    parser.add_argument("--out_dir", dest='out_dir', default="output", type=str, help="Directory, which output files are written to")
    parser.add_argument("--k", dest='k', default="1", type=int, help="k of NN")
    opts = parser.parse_args()
    return opts


if __name__ == '__main__':

    start = time.time()

    # Read user input
    opts = readUserArguments()

    target_files = glob.glob(opts.target)

    out_dir = opts.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    source_fname, source_ext = os.path.splitext(opts.source)

    allowed_extensions = [".las", ".laz"]

    # open source cloud (check if las or ascii) and get header text of each cloud
    if source_ext == ".las" or source_ext == ".laz":
        inFile = File(opts.source, mode="r")
        source_cloud = readLas(inFile)
    else:
        print("File format not supported. Please use target files with extensions '.las' or '.laz'.")
        sys.exit()

    # get number of source points
    n_source = source_cloud.shape[0]
    print("%i points in source cloud read." % n_source)

    # Take only XYZ of source cloud for 3D KDTree
    kd_tree = KDTree(source_cloud[:, :3], leafsize=64)

    if opts.has_header == 0:
        header = None
        header_txt = None

# keep points from source cloud (ALS/ULS) only if there is a point in target (TLS) within a certain distance
    for index, target in enumerate(target_files):

        fname_target, target_ext = os.path.splitext(os.path.basename(target))
        
        # open target cloud
        if target_ext == ".las" or target_ext == ".laz":
            inFile_target = File(target, mode="r")
            target_cloud = readLas(inFile_target)
        else:
            print("File format not supported. Please use target files with extensions '.las' or '.laz'.")
            sys.exit()
        
        # Get number of points
        n_target = target_cloud.shape[0]
        print("%i points in target cloud %i read." % (n_target, index+1))

        for radius in [0.1, 0.3, 0.5, 1.0]:
            # dist - stores the Euclidean linear distances to all neighbours
            # idx  - stores the array index of the neigbhour in  the original data array: needed to access the point data itself
            dist, idx = kd_tree.query(target_cloud[:, :3], k=opts.k, sqr_dists=False, eps=0.0, distance_upper_bound=radius)
            
            idx = np.unique(idx)
            idx = idx[idx != n_source]
            
            out_dir = os.path.join(opts.out_dir, "search_radius_" + str(int(radius*100)))
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            if source_ext == ".laz":
                output_las = out_dir + "\\" + fname_target + "_ULS.laz"
            else:
                output_las = out_dir + "\\" + fname_target + "_ULS.las"

            # write result to las file
            print("Write result to %s" % output_las)
            outfile = File(output_las, mode="w", header=inFile.header)
            outfile.points = inFile.points[idx]
            outfile.close()

    end = time.time()
    runtime = end - start
    minutes = math.floor((end - start) / 60)
    seconds = runtime - 60 * minutes
    print("Runtime: %i minutes and %.3fs seconds" % (minutes, seconds))
