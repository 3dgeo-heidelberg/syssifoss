#!/bin/python
# -*- coding: utf-8 -*-

"""
least squares ellipse fitting using code from https://github.com/ndvanforeest/fit_ellipse by user ndvanforeest
The DBH (as the mean of the two axes of the ellipse) is written to the file and the stem points with the
fitted ellipse is plotted and saved to a file (output files named by tree_ID taken from the filename)

Usage: python DBH_ellipse.py <pointcloud.las> <tree_position.txt> <out_dir>
    infile.las          is a LAS/LAZ file of the tree (multiple infiles can be specified using wildcards)
    tree_position.txt   is an ASCII file with the tree position (comma-separated x,y,z)
                        - the same tree ID has to be contained in the filename of infile.las and translation.txt
                        - the tree ID is expected in the format: speciesID_plotID_number
                        - the z-value from tree position.txt is used to normalize the height of the point cloud
    out_dir             is the directory to write the output files to. One txt-file is generated per tree

wildcards can be used to select multiple files (e.g. data\*.laz)

Example usage:
"""


from laspy.file import File
import sys
import glob
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from fit_ellipse import *


in_file = sys.argv[1]
translation = sys.argv[2]
out_dir = sys.argv[3]

infiles = glob.glob(in_file)
translation_files = glob.glob(translation)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

for i, fid in enumerate(infiles):

    if i % 10 == 0:
        if i == 0:
            print("Start processing...")
        else:
            print("Processing point cloud %d..." % i)

    filebase = os.path.basename(fid)
    filename = os.path.splitext(filebase)[0]
    fileparts = filename.split("_")
    tree_id = "%s_%s_%s" % (fileparts[0], fileparts[1], fileparts[2])
    outf = os.path.join(out_dir, "%s_m_%s_%s" % (tree_id, fileparts[3], fileparts[5]) + "DBH_ellipse")
    print(outf)

    # load point cloud
    inFile = File(fid, mode="r")
    pc = np.array([inFile.x, inFile.y, inFile.z]).T
    inFile.close()

    # load translation - has to contain the tree-id of the input file!
    t_file = [s for s in translation_files if tree_id in s]
    with open(t_file[0], "r") as fileobj:
        t = float(fileobj.read().split(",")[2])
    pc[:, 2] -= t

    # define breast height and height of stem slice
    breast_height = 1.3
    slice_height = 0.04

    # extract stem slice
    stem_slice = pc[(pc[:, 2] >= breast_height - slice_height / 2) & (pc[:, 2] <= breast_height + slice_height / 2)]

    #prepare plot
    fig, ax = plt.subplots(figsize=(10,10))
    ax.scatter(stem_slice[:, 0], stem_slice[:, 1])
    plt.axis("equal")

    x = stem_slice[:, 0]
    y = stem_slice[:, 1]
    x_mean = x.mean()
    y_mean = y.mean()
    x = x - x_mean
    y = y - y_mean
    a, b, center0, center1, phi = fit_ellipse(x, y)
    center0 += x_mean
    center1 += y_mean
    center = (center0, center1)
    a = a*2
    b = b*2
    ell = Ellipse(center, a, b, phi * 180 / np.pi, color="r", lw = 4.0, fill=False)  # convert phi from degree to rad

    # diameter as mean of the major and minor axes of the ellipse
    diameter = (a+b)/2*100 # convert from m to cm

    # write some info to the title
    plt.title("%s; slice: %.2f-%.2f m; a = %.2f, b = %.2f"
              % (filename[:9], breast_height-slice_height/2, breast_height+slice_height/2,  a, b))
    ax.add_artist(ell)
    plt.savefig(outf+".png", dpi=100)
    plt.clf()
    plt.close()

    # write DBH to output file
    with open(outf+".txt", "w") as outfile_txt:
        outfile_txt.write("DBH_cm\n%.1f\n" % diameter)
