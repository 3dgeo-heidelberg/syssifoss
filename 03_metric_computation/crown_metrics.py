#!/bin/python
# -*- coding: utf-8 -*-
# Hannah Weiser (Heidelberg University)
# h.weiser@stud.uni-heidelberg.de

"""
Driver for the computation of tree height, crown base height, crown area, diameter at breast height

Usage: python tree_metrics.py <infile.las> <tree_position.txt> <out_dir> <concave_exe> (<metric_file.txt>)
    infile.las          is a LAS/LAZ file of the tree
    tree_position.txt   is an ASCII file with the tree position (comma-separated x,y,z)
                        - the same tree ID has to be contained in the filename of infile.las and translation.txt
                        - the tree ID is expected in the format: speciesID_plotID_number
                        - the z-value from tree position.txt is used to normalize the height of the point cloud
    out_dir             is the directory to write the output files to. One txt-file is generated per tree
    metric.txt          is a metrics text file, containing a DBH value in the column "DBH_cm" (OPTIONAL)

wildcards can be used to select multiple files (e.g. data\*.laz)

Example usage:

"""

import numpy as np
import pandas as pd
import TreeMetrics
from laspy.file import File     # laspy has to be installed: pip install laspy
import sys
import matplotlib.pyplot as plt
import os
import glob
from scipy.spatial import ConvexHull, convex_hull_plot_2d


def height_percentiles(point_cloud, qq=[]):
    """
    @brief Function to compute height quantiles for a point cloud.

    @param point_cloud: numpy array of point cloud, e.g. (n, 3) point cloud with n = number of points;
                        z-coordinate is expected at index 2
    @param qq:  List of percent ranks to compute quantiles for, e.g. [95] for the 95th percentile
                If this parameter is not provided, the deciles will be computed

    @return: List of quantiles
    """

    if qq != []:
        pass
    else:
        qq = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90])

    h_q = np.zeros(len(qq), dtype=float)
    for idx, qq in enumerate(qq):
        h_q[idx] = np.percentile(point_cloud[:, 2], qq)

    return h_q


in_file = sys.argv[1]
translation = sys.argv[2]
out_dir = sys.argv[3]
# if DBH measurements available, load them from here:
if len(sys.argv) > 4:
    metrics_files = glob.glob(sys.argv[4])
else:
    metric_files = None
concave_exe = "concave.exe"

infiles = glob.glob(in_file)
translation_files = glob.glob(translation)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

for fid in infiles:

    filebase = os.path.basename(fid)
    filename = os.path.splitext(filebase)[0]
    fileparts = filename.split("_")
    tree_id = "%s_%s_%s" % (fileparts[0], fileparts[1], fileparts[2])

    # load point cloud
    inFile = File(fid, mode="r")
    pc = np.array([inFile.x, inFile.y, inFile.z])
    pc = pc.T

    # load translation - has to contain the tree-id of the respective infile
    t_file = [s for s in translation_files if tree_id in s]
    with open(t_file[0], "r") as fileobj:
        t = float(fileobj.read().split(",")[2])
    # normalize z-values by the height of the tree position
    pc[:, 2] -= t

    print("Finished reading input cloud %s" % filebase)

    # if DBH is available, load it with this code:
    if metrics_files is not None:
        metric_file = [f for f in metrics_files if tree_id in f]
        if not metric_file:
            DBH = 0.0
            continue
        df = pd.read_csv(metric_file[0], sep=" ")
        DBH = df.loc[0, "DBH_cm"]
    else:
        DBH = 0.0

    out_path = os.path.join(out_dir, "%s_m_%s_%s.txt" % (tree_id, fileparts[3], fileparts[5]))

    # open attribute output file
    with open(out_path, "w") as outfile:
        # if no DBH was available from other measurement and the tree is taller than 15 m, set it to 0.5
        if np.max(pc[:, 2]) > 15 and DBH == 0.0:
            DBH = 0.5

        # write header of output file
        outfile.write("height_m "
                      "crown_base_height_m "
                      "crown_projection_area_convex_hull_m2 "
                      "crown_projection_area_concave_hull_m2 "
                      "mean_crown_diameter_m\n")

        tm = TreeMetrics.TreeMetrics(pc)

        print("Tree height is calculated...")
        max_h = np.max(pc[:, 2])
        h99_9 = height_percentiles(pc, [99.9])[0]
        # use the 99.9 percentile as height if the 99.9 percentile is more than 1.0 m below h_max (due to an outlier)
        # else use max_h
        if (max_h - h99_9) > 1.0:
            h = h99_9
        else:
            h = max_h

        print("compute_cbh is computed...")
        branch_length = 1.0
        lim, max_dists, CBH = tm.compute_cbh(pc, 0.1, branch_length, DBH)

        # plot the compute_cbh
        plt.plot(max_dists, lim, '-g')
        plt.axhline(y=CBH, color='y', linestyle='-', linewidth=1.)
        plt.xlabel("Max. horizontal point distance [m]")
        plt.ylabel("Center of height bin [m]")
        plt.axis('square')
        figure = plt.gcf()
        figure.set_size_inches(4, 4)
        plt.savefig(out_path.replace(".txt", "_CBH.png"))
        plt.clf()
        plt.close('all')
        print("compute_cbh plot is saved to: %s" % out_path.replace(".txt", "_CBH.png"))

        print("Crown area is computed with convex hull algorithm...")

        # only use points above the compute_cbh
        pc_crown = pc[pc[:, 2] >= CBH]

        # except if then there are less than 3 points left
        if pc_crown.shape[0] <= 3:
            pc_crown = pc
        Hull = ConvexHull(pc_crown[:, :2])
        convex_hull_poly = pc_crown[:, :2][Hull.vertices]
        CA_convex_hull = tm._polygon_area(convex_hull_poly)

        # plot the convex hull
        convex_plt = convex_hull_plot_2d(Hull)
        plt.axis("image")
        x_ax_min, x_ax_max, y_ax_min, y_ax_max = plt.axis()
        plt.title("Convex Hull")
        convex_plt.set_size_inches(5., 5.)
        convex_plt.savefig(out_path.replace(".txt", "_convex_hull.png"), dpi=90)
        convex_plt.clf()
        plt.close('all')

        print("Crown area is computed with concave hull algorithm...")

        # save point cloud as txt-file in a temporary folder to run concave hull algorithm
        temp_dir = out_dir + "\\temp\\"
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        np.savetxt(temp_dir + filename + ".txt", pc_crown[:, :2])
        CA_concave_hull, concave_hull_poly = tm.ca_concave(concave_exe, temp_dir + filename + ".txt")

        if concave_hull_poly.shape[0] != 0:
            # plot the concave hull
            plt.ylim(y_ax_min, y_ax_max)
            plt.xlim(x_ax_min, x_ax_max)
            concave_plt = tm.plot_hull(pc_crown, concave_hull_poly)
            plt.title("Concave hull")
            plt.axis("image")
            concave_plt.savefig(out_path.replace(".txt", "_concave_hull.png"), dpi=90)
            concave_plt.clf()

            # if a concave hull was successfully computed, derive the crown diameter from the concave hull
            print("Crown diameter is computed...")
            max_crown_width, max_cross_width, average_crown_width, crown_width_plt = tm.crown_diameter(concave_hull_poly)
            crown_width_plt.savefig(out_path.replace(".txt", "_crown_width.png"), dpi=90)
            crown_width_plt.clf()
            plt.close('all')
        else:
            # use the convex hull instead for deriving the crown diameter
            print("Crown diameter is computed...")
            max_crown_width, max_cross_width, average_crown_width, crown_width_plt = tm.crown_diameter(convex_hull_poly)
            crown_width_plt.savefig(out_path.replace(".txt", "_crown_width.png"), dpi=90)
            crown_width_plt.clf()
            plt.close('all')

        # write everything to the output file
        outfile.write("%.3f %.3f %.3f %.3f %.3f" % (h, CBH, CA_convex_hull, CA_concave_hull, average_crown_width))

