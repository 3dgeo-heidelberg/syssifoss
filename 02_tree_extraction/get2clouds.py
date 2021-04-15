#!/bin/python
# -*- coding: utf-8 -*-

"""
Grab points and attributes from source point cloud using points from target point cloud:
  - For each point in the target point cloud, the closest k points in the source point cloud within a specified
    search_radius are written to the output point cloud (with attributes from source point cloud)
  - Assignment is done via Euclidean distance (provided by the search_radius parameter)

Example use case: Extract single tree point cloud from a TLS plot point cloud with all attributes using
(downsampled) segmented TLS point clouds of single trees (only xyz)

Important:
- One source point cloud and multiple target point clouds can be provided.
- Source and target point clouds can be in LAS/LAZ format or ASCII format. The output point cloud will be in the same
  format as the source point cloud

(c) H Weiser, 2020 from Bernhard Hoefle, 2019
"""

import sys
import argparse
from pykdtree.kdtree import KDTree  # Install: https://github.com/storpipfugl/pykdtree
import numpy as np
import pandas as pd  # pandas needs to be installed: pip install pandas
import time
from laspy.file import File  # laspy needs to be installed: pip install laspy
import math
import os
import glob


def getHeaderLine(filename):
    with open(filename, 'r') as fp:
        header_txt = fp.readline().strip()
    return header_txt


def getDimensions(file):
    dimensions = ""
    for dim in file.point_format:
        dimensions += " " + dim.name
    return dimensions


def readLas(file):
    dimensions = getDimensions(file)
    source_cloud = np.array([file.x, file.y, file.z]).T
    return source_cloud, dimensions


def readUserArguments():
    # Generate help and user argument parsing
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="(c) Bernhard Hoefle (2019) - 3DGeo @ Heidelberg University")
    parser.add_argument("--source", dest='source', type=str,
                        help="Source point cloud from which attributes and points shall be taken. It can be in ascii or las format",
                        required=True)
    parser.add_argument("--target", dest='target', type=str,
                        help="Target point cloud(s) to which the attributes from source point cloud should be appended and which is the spatial select: only points of source close to points of the target are considered",
                        required=True)
    parser.add_argument("--out_dir", dest='out_dir', default="output", type=str,
                        help="Directory, to which output files are written")
    parser.add_argument("--k", dest='k', default="1", type=int, help="k of NN")
    parser.add_argument("--search_radius", dest='search_radius', type=float, default=0.1,
                        help="Maximum 3D Euclidean distance to take point from source")
    parser.add_argument("--sep", dest='sep', type=str, default=" ", help="field/column separator")
    parser.add_argument("--has_header", dest='has_header', type=int, default="1",
                        help="Do the input files contain a header line? (1 ... yes / 0 ... no / -1 yes but do not use = skip first line)")
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

    allowed_extensions = [".las", ".laz", ".asc", ".xyz", ".txt"]

    # open source cloud (check if las or ascii) and get header text of each cloud
    if source_ext not in allowed_extensions:
        print(
            "File format not supported. Please use target files with extensions '.las', '.laz', '.asc', '.xyz' or '.txt'")
        sys.exit()
    elif source_ext == ".las" or source_ext == ".laz":
        inFile = File(opts.source, mode="r")
        source_cloud, header_txt_source = readLas(inFile)
    elif opts.has_header == -1:
        source_cloud = pd.read_csv(opts.source, opts.sep, header=0, skipinitialspace=True, skiprows=0).values
    else:
        source_cloud = pd.read_csv(opts.source, opts.sep, header=opts.has_header, skipinitialspace=True).values
        if opts.has_header == 1:
            header_txt_source = getHeaderLine(opts.source)

    # get number of source points
    n_source = source_cloud.shape[0]
    print("%i points in source cloud read." % (n_source))

    # Take only XYZ of source cloud for 3D KDTree
    kd_tree = KDTree(source_cloud[:, :3], leafsize=64)

    if opts.has_header == 0:
        header = None
        header_txt = None

    # keep points from source cloud only if there is a point in target within a certain distance
    for index, target in enumerate(target_files):

        fname_target, target_ext = os.path.splitext(os.path.basename(target))

        # Handle header for pandas reading (0... take first line as header; None.. no header line in input)
        if opts.has_header == 1 and (target_ext != ".las" and target_ext != ".laz"):
            header = 0
            header_txt_target = getHeaderLine(target)

        # open target cloud
        if target_ext == ".las" or target_ext == ".laz":
            inFile_target = File(target, mode="r")
            target_cloud, header_txt_target = readLas(inFile_target)
        else:
            target_cloud = pd.read_csv(target, opts.sep, header=header).values

        # Get number of points
        n_target = target_cloud.shape[0]
        print("%i points in target cloud %i read." % (n_target, index + 1))

        # dist - stores the Euclidean linear distances to all neighbours
        # idx  - stores the array index of the neigbhour in  the original data array: needed to access the point data itself
        dist, idx = kd_tree.query(target_cloud[:, :3], k=opts.k, sqr_dists=False, eps=0.0,
                                  distance_upper_bound=opts.search_radius)

        idx = np.unique(idx)
        idx = idx[idx != n_source]

        if source_ext == ".asc" or source_ext == ".txt" or source_ext == ".xyz" or source_ext == ".csv":
            output_ascii = opts.out_dir + "\\" + fname_target + ".xyz"
            try:
                # respective points from source point cloud
                output_cloud = source_cloud[idx]
                print("Write result to %s" % output_ascii)
                # Write result which is VERY SLOW like this
                if opts.has_header == 1:
                    np.savetxt(output_ascii, output_cloud, fmt="%.3f", comments='', delimiter=opts.sep,
                               header=header_txt_source)
                elif opts.has_header == 0 or opts.has_header == -1:
                    np.savetxt(output_ascii, output_cloud, fmt="%.3f", comments='', delimiter=opts.sep)

            except Exception as err:
                print("No corresponding points were found. Consider increasing search_radius")
                print("Python Error: ", err)

        elif source_ext == ".las" or source_ext == ".laz":
            if target_ext == ".las":
                output_las = opts.out_dir + "\\" + fname_target + ".las"
            elif target_ext == ".laz":
                output_las = opts.out_dir + "\\" + fname_target + ".laz"
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
