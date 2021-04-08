#!/bin/python
# -*- coding: utf-8 -*-
# Hannah Weiser, Heidelberg University, July 2020 (code adapted from:
# course: machine learning for geographic applications 2020
# exercise 7: robust methods, Lukas Winiwarter)

"""
robust circle fitting using RANSAC - works even for noisy point clouds (e.g. stem slices with additional leaf points)
The DBH is written to the file and the stem points with the fitted circle is plotted
and saved to a file (output files named by tree_ID taken from the filename, can be changed in line 68)

Usage: python DBH_RANSAC.py <pointcloud.las> <tree_position.txt> <out_dir>
    infile.las          is a LAS/LAZ file of the tree (multiple infiles can be specified using wildcards)
    tree_position.txt   is an ASCII file with the tree position (comma-separated x,y,z)
                        - the same tree ID has to be contained in the filename of infile.las and translation.txt
                        - the tree ID is expected in the format: speciesID_plotID_number
                        - the z-value from tree position.txt is used to normalize the height of the point cloud
    out_dir             is the directory to write the output files to. One txt-file is generated per tree

wildcards can be used to select multiple files (e.g. data\*.laz)

Example usage:
"""

from laspy.file import File # pip install laspy
import numpy as np
import sys
import glob
import os
import matplotlib.pyplot as plt

def is_in_circle(x, y , xm, ym, r, eps=0.01):
    rad2 = (x-xm)**2 + (y-ym)**2
    cond1 = rad2 <= (r+eps)**2
    cond2 = rad2 >= (r-eps)**2
    return np.logical_and(cond1, cond2)

def circle_from_points(p1,p2,p3):
    p1p2 = p2-p1
    p2p3 = p3-p2
    cent12 = (p1+p2)/2
    cent23 = (p2+p3)/2
    n12 = np.array([p1p2[1], -p1p2[0]])
    n23 = np.array([p2p3[1], -p2p3[0]])
    b = cent23-cent12
    A = np.array([n12, n23])
    x = np.linalg.pinv(A).dot(b)
    C = x[0] * n12 + cent12
    r = np.linalg.norm(p1-C)
    return C[0], C[1], r


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

    filebase = os.path.basename(fid)
    filename = os.path.splitext(filebase)[0]
    outf = os.path.join(out_dir, "%s_m_%s_%s" % (tree_id, fileparts[3], fileparts[5]) + "DBH_RANSAC")
    print(outf)

    # load point cloud
    inFile = File(fid, mode="r")
    pc = np.array([inFile.x, inFile.y, inFile.z])
    pc = pc.T
    inFile.close()

    # load translation; has to contain the tree-id of the input file!
    t_file = [s for s in translation_files if tree_id in s]
    with open(t_file[0], "r") as fileobj:
        t = float(fileobj.read().split(",")[2])
    pc[:, 2] -= t

    print("Finished reading input cloud %s" % filebase)
    breast_height = 1.3
    slice_height = 0.04

    stem_slice = pc[(pc[:, 2] >= breast_height - slice_height / 2) & (pc[:, 2] <= breast_height + slice_height / 2)]
    if stem_slice.shape[0] > 2:
        fig, ax = plt.subplots(figsize=(10,10))
        ax.scatter(stem_slice[:, 0], stem_slice[:, 1])
        plt.axis("equal")
        x = stem_slice[:, 0]
        y = stem_slice[:, 1]

        x_copy = x
        y_copy = y

        final_score = 0
        color = 'r'
        score_max = 0
        circ_max = None
        i = 0
        while i < 1000:
            try:
                three_points = np.random.choice(range(len(x_copy)), 3, replace=False)
            except ValueError:
                break
            x3 = x_copy[three_points]
            y3 = y_copy[three_points]
            p1 = np.array([x3[0], y3[0]])
            p2 = np.array([x3[1], y3[1]])
            p3 = np.array([x3[2], y3[2]])
            cx, cy, r = circle_from_points(p1, p2, p3)
            if r > 0.8:
                continue
            score = np.count_nonzero(is_in_circle(x_copy, y_copy, cx, cy, r))
            if score > score_max:
                score_max = score
                circ_max = cx, cy, r
                n_max = x_copy.shape[0]
                circle1 = plt.Circle((circ_max[0], circ_max[1]), circ_max[2], color="yellow", lw=1.0, fill=False)
                ax.add_artist(circle1)
            i += 1

        # plot circle
        circle1 = plt.Circle((circ_max[0], circ_max[1]), circ_max[2], color=color, lw = 4.0, fill=False)
        diameter = (circ_max[2]*2)*100
        ax.add_artist(circle1)
        plt.title("%s; slice: %.2f-%.2f m; Final score: %i/%i; r = %.3f m" % (filename[:9], breast_height-slice_height/2, breast_height+slice_height/2, score_max, n_max, circ_max[2]))
        plt.savefig(outf+".png", dpi=100)
        plt.clf()
        plt.close()

    else:
        diameter = np.nan

    with open(outf+".txt", "w") as outfile_txt:
        outfile_txt.write("DBH_cm\n%.1f\n" % diameter)
