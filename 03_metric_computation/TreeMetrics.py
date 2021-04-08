#!/bin/python
# -*- coding: utf-8 -*-
# Hannah Weiser (Heidelberg University)
# h.weiser@stud.uni-heidelberg.de

"""
Class for computation of standard tree metrics

point_cloud        (n,3) np.array with columns [x, y, nZ]
"""

import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist, squareform
import math
import matplotlib.pyplot as plt
import subprocess
import warnings
import os

class TreeMetrics():

    def __init__(self, point_cloud):
        self.pc = point_cloud

    def ca_concave(self, exe, file):
        """
        @brief Function to compute the concave hull:

        @param exe: path of the C++ program using the concave hull algorithm by Moreira & Santos (2007), which uses
                    a k nearest neighbours approach
                    downloaded: https://www.codeproject.com/Articles/1201438/The-Concave-Hull-of-a-Set-of-Points
        @param file: ASCII-file with the point cloud;
                    first two fields should be x and y;
                    fields can be comma, tab or space delimited

        @return: [area of the concave hull, Ordered sequence of the vertices of the hull polygon]
        """

        hull_file = file.replace(".txt", "_hull.xyz")
        # if the hull file already exists (i.e. has been previously computed), do not compute again
        if not os.path.exists(hull_file):
            subprocess.call([exe, file, "-out", hull_file])
        hull_polygon = np.loadtxt(hull_file)
        if hull_polygon.shape[0] != 0:
            ca_concave_hull = self._polygon_area(hull_polygon)
        else:
            ca_concave_hull = np.nan

        return ca_concave_hull, hull_polygon

    def plot_hull(self, pc, hull_pts):
        """
        @brief Function to plot the (concave or convex) hull

        @param pc: numpy array of point cloud for which the hull was computed, e.g., of shape (n, 3)
        @param hull_pts: Ordered sequence of the vertices of the hull polygon
        @return: the matplotlib plot
        """

        plt.close("all")

        x = pc[:, 0]
        y = pc[:, 1]

        hull_x = hull_pts[:, 0]
        hull_y = hull_pts[:, 1]

        fig = plt.figure(1, figsize=(5, 5), dpi=90)
        ax = fig.add_subplot(111)
        ax.plot(x, y, 'o', color="green")
        ax.plot(hull_x, hull_y, color="#ed7d31", lw=3)
        plt.axis("image")

        return plt

    def crown_diameter(self, hull_poly):
        """
        @brief Function to compute the crown diameter from a hull polygon
        The diameter is defined as the mean of the maximum crown width and the maximum cross width

        @param hull_poly: Ordered sequence of the vertices of the hull polygon

        @return:    [The maximum crown width
                    The maximum cross width
                    The (mean) crown diameter
                    The plot]
        """

        # sample regularly spaced points on the dripline of the tree (= concave or convex hull)
        hull_poly = np.array(hull_poly)
        points = self._sample_points_on_polygon_outline(hull_poly, spacing=0.05)

        # distance matrix
        d = pdist(points)
        d = squareform(d)

        # get maximum crown width and corresponding points
        max_crown_width = np.nanmax(d)
        i_a, i_b = np.unravel_index(np.argmax(d), d.shape)
        a, b = points[i_a], points[i_b]
        # vector between the points furthest away from each other
        max_vec = b - a

        # iterate through the distance matrix, get the vector between respective points and check if it is perpendicular
        # to max_vec (two vectors are perpendicular if their dot product is zero). Update the found points and the
        # curr_highest_dist if it is higher than the ones previously found
        curr_highest_dist = 0
        for i, dist in np.ndenumerate(d):
            curr_p1, curr_p2 = points[i[0]], points[i[1]]
            vec = curr_p2 - curr_p1
            dot = np.sum(vec*max_vec)
            if round(dot) == 0 and dist > curr_highest_dist:
                curr_highest_dist = dist
                p1, p2 = curr_p1, curr_p2

        # Plotting
        # plot a line for the maximum width and the maximum cross width
        plt.plot([a[0], b[0]], [a[1], b[1]], marker="*", color="g")
        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], marker="*", color="r")
        # plot the points defining the dripline (the hull)
        plt.scatter(points[:, 0], points[:, 1])
        plt.axis("image")

        return max_crown_width, curr_highest_dist, (max_crown_width+curr_highest_dist)/2, plt

    def cv_convex_hull(self, cbh):
        """
        @brief Function to compute the crown volume using the convex hull algorithm

        @param cbh: the height of the crown base [m]

        @return: the crown volume [m2]
        """

        # get all points above the cbh (= crown points)
        crown_pc = self.pc[self.pc[:,2] >= cbh]
        # compute the convex hull and calculate its volume
        crown_volume = ConvexHull(crown_pc).volume

        return crown_volume

    def cv_concave_hull_adaptive(self, t0, CBH, exe, temp_fn):
        # calculates the crown volume using concave hull areas of crown slices with adaptive thickness (Yan et al. 2019)
        # exe is the path of the C++ program to the concave hull algorithm using a k nearest neighbours approach
        # by Moreira & Santos (2007)
        # (can be downloaded here: https://www.codeproject.com/Articles/1201438/The-Concave-Hull-of-a-Set-of-Points)

        warnings.filterwarnings('error')
        np.errstate(all='raise')

        # slice crown into slices with uniform thickness t0
        nZ = self.pc[:, 2]
        lim = np.arange(CBH, np.ceil(np.max(nZ)*10)/10 + t0, t0)

        # catch slice with less than three points and merge with upper adjacent slice until there are at least 3 points
        pCount = np.histogram(nZ, lim)[0]
        while min(pCount) < 3:
            for i, count in enumerate(pCount):
                if count < 3:
                    lim = np.delete(lim, i+1)
                    pCount = np.histogram(nZ, lim)[0]
                    break

        # get concave area of each slice
        areas = np.zeros(pCount.shape)
        for i, limit in enumerate(lim):
            if i == lim.shape[0]-1:
                break
            pc_slice = self.pc[(self.pc[:, 2] >= limit) & (self.pc[:,2] < lim[i+1])]
            np.savetxt(temp_fn, pc_slice)
            hull_file = temp_fn.replace(".xyz", "_hull.xyz")
            subprocess.call([exe, temp_fn, "-out", hull_file])
            if os.stat(hull_file).st_size == 0:
                try:
                    # if no concave hull can be found (can be the case for the models due to the voxel structure), calculate
                    # convex hull for this section instead
                    hull_polygon = pc_slice[ConvexHull(pc_slice).vertices]
                    CA_concave_hull = self._polygon_area(hull_polygon)
                except: # exception will be raised, when less than 3 points in slice - should not happen but idk..
                    CA_concave_hull = 0
            else:
                hull_polygon = np.loadtxt(hull_file)
                CA_concave_hull = self._polygon_area(hull_polygon)
            areas[i] = CA_concave_hull

        Smean = np.mean(areas)
        std = np.std(areas)
        # categorize according to area using the formula
        cat = np.zeros(pCount.shape)
        for i, S in enumerate(areas):
            if S - Smean >= 0:
                C = ((S - Smean)/std) + 1
            elif S - Smean < 0:
                C = ((S - Smean) / std) - 1
            cat[i] = C
        cat = np.trunc(cat)

        # adapt limits, so that sections of the same category are merged
        lim_adaptive = np.zeros(lim.shape)
        curr_slice = 1
        curr_cat = cat[0]
        lim_adaptive[0] = lim[0]
        for i, a in enumerate(cat):
            if a != curr_cat:
                curr_slice += 1
                curr_cat = a
                lim_adaptive[curr_slice-1] = lim[i]
            if i == cat.shape[0]-1:
                lim_adaptive[curr_slice] = lim[-1]
                print("done")
        lim_adaptive = lim_adaptive[:curr_slice+1]

        # compute concave hull for each section and get slice height
        areas_adapted = np.zeros(curr_slice+1)
        heights = np.zeros(curr_slice+1)
        for i, limit in enumerate(lim_adaptive):
            if i == curr_slice:
                break
            pc_slice = self.pc[(self.pc[:, 2] >= limit) & (self.pc[:, 2] < lim_adaptive[i+1])]
            slice_height = lim_adaptive[i+1]-limit
            np.savetxt(temp_fn, pc_slice)
            hull_file = temp_fn.replace(".xyz", "_hull.xyz")
            subprocess.call([exe, temp_fn, "-out", hull_file])
            try:
                hull_polygon = np.loadtxt(hull_file)
            except UserWarning:
                # if no concave hull can be found (can be the case for the models due to the voxel structure), calculate
                # convex hull for this section instead
                hull_polygon = pc_slice[ConvexHull(pc_slice).vertices]
            CA_concave_hull = self._polygon_area(hull_polygon)
            areas_adapted[i] = CA_concave_hull
            heights[i] = slice_height
        print(heights)
        # compute volume as the mean of the area of the two planes limiting a section... (not sure about this step)
        volumes = np.zeros(curr_slice)
        for i in range(curr_slice):
            volumes[i] = (areas_adapted[i]+areas_adapted[i+1])/2 * heights[i]

        # sum of all section volumes
        return np.sum(volumes)

    def compute_cbh(self, pts, dZ, threshold, dbh):
        """
        @brief Function to compute the crown base height (CBH) by splitting the point cloud into height bins and
        deriving the minimum height in which the maximum x-y distance between any two points exceeds a certain threshold

        @param pts: tree point cloud
        @param dZ: size of height bins to split the point cloud into
        @param threshold:   Rough branch size to consider large enough to count as the crown base;
                            x-y distance in height bin has to exceed DBH + threshold
        @param dbh: the DBH of the tree

        @return:    [The height bin centers,
                    The maximum x-y distances in each height bin,
                    The CBH]

        """
        # split point cloud in bins of certain height
        lim, bin_centers = self._split_into_bins(pts, dZ)
        cbh = 0.0

        # maximum point distances for each bin
        max_dists = np.zeros(len(bin_centers))
        for i, limit in enumerate(lim):
            if i != 0:
                pc = pts[(pts[:, 2] <= limit) & (pts[:, 2] > lim[i - 1])]
            else:
                pc = pts[pts[:, 2] <= lim[i+1]]
            if pc.shape[0] > 2:
                hull = ConvexHull(pc[:, :2])
                d = pdist(pc[:, :2][hull.vertices])
                d = squareform(d)
                dist = np.nanmax(d)
                max_dists[i] = dist

        # starting from the bottom, get the first height bin, where the maximum point distance exceeds DBH + threshold
        for i, dist in enumerate(max_dists):
            start = 0.0
            # first check to find first height bin which does not belong ground points (if there are any)
            # this step is only necessary, if the point cloud contains ground points
            if dist <= 1.3 and i > 10: # 1.3, 10 ALS: 1.3, 5
                start = i
                break
        # start = 0.0
        for i, dist in enumerate(max_dists):
            # check if threshold is exceeded and if distance is higher than distance in previous height bin
            if i > start and dist >= (threshold + round(dbh, 2)) and dist > max_dists[i - 1]:
                cbh = bin_centers[i]
                break

        return bin_centers, max_dists, cbh

    def CBH_supervised(self, pts, dZ):
        """
        @brief: Function to detect the CBH interactively -> the user clicks to the height, where he/she assumes the CBH
                and the respective z value is stored as CBH

        @param pts: tree point cloud
        @param dZ: size of height bins to split the point cloud into

        @return:    [The height bin centers,
                    The maximum x-y distances in each height bin,
                    The CBH]
        """

        # split point cloud in bins of certain height
        lim, bin_centers = self._split_into_bins(pts, dZ)

        # get the maximum point distance of each height bin
        max_dists = np.zeros(len(bin_centers))
        for i, limit in enumerate(lim):
            if i == 0:
                pc = pts[pts[:, 2] <= lim[i+1]]
            else:
                pc = pts[(pts[:, 2] <= limit) & (pts[:, 2] > lim[i-1])]
            if pc.shape[0] > 2:
                d = pdist(pc[:, :2])
                d = squareform(d)
                max_dist = np.nanmax(d)
                max_dists[i] = max_dist

        plt.plot(max_dists, bin_centers, '-g')
        plt.xlabel("Max. horizontal point distance [m]")
        plt.ylabel("Center of height bin [m]")
        plt.axis('square')
        print("Select the CBH")

        # here, the user does his job
        cbh = plt.ginput(1, timeout=-1)[0][1]
        plt.close()
        plt.clf()

        return bin_centers, max_dists, cbh

    # ----------------------------------------------------------------------
    # class methods

    def _polygon_area(self, points):
        """
        @brief  Shoelace formula
                Returns the area of the polygon whose vertices are given (ordered) by the
                sequence points.

        @param points: Ordered sequence of the vertices of the polygon

        @return: The area of the polygon
        """

        area = 0
        q = points[-1]
        for p in points:
            area += p[0] * q[1] - p[1] * q[0]
            q = p

        return abs(area / 2)


    def _sample_points_on_polygon_outline(self, poly, spacing=0.05):
        """
        samples evenly spaced points on the outline of a polygon

        :param poly: numpy array containing the points defining a polygon (one point per row)
        :param spacing: number to control the point spacing
        :return: numpy array containing evenly spaced points on the polygon outline
        """
        for i in range(poly.shape[0]):
            vec = poly[i] - poly[i - 1]
            dist = math.sqrt(vec[0] ** 2 + vec[1] ** 2)
            parts = round(dist / spacing)
            if parts != 0:
                new_points = list(zip(np.linspace(poly[i - 1][0], poly[i][0], parts),
                                      np.linspace(poly[i - 1][1], poly[i][1], parts)))
            else:
                new_points = poly[i]
            if i == 0:
                points = np.array(new_points)
            else:
                new_points = np.array(new_points)
                points = np.vstack((points, new_points))

        return points


    def _split_into_bins(self, pts, dZ):
        """
        Splits a point cloud into bins of a certain height dZ and return the bin limits and center
        :param pts: 3D point cloud (shape (n, 3) with n = number of points and columns x, y and normalized z)
        :param dZ: bin height
        :return: limits, centers
        """

        # get (normalized) height
        nZ = pts[:, 2]

        # define limits of bins within which points are counted
        num = math.ceil(np.max(nZ) / dZ)
        upper_lim = dZ * num
        lim = np.arange(0, upper_lim, dZ)
        # define bin centers
        bin_centers = np.arange(0 + dZ / 2, upper_lim + dZ / 2, dZ)

        return lim, bin_centers
