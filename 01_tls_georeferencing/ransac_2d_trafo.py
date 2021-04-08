#!/bin/python
# -*- coding: utf-8 -*-
__author__ = "Lukas Winiwarter"
__copyright__ = "Copyright 2019, 3DGeo/Uni Heidelberg"
__license__ = "MIT"
__version__ = "1.0"
__email__ = "lukas.winiwarter@uni-heidelberg.de"
__status__ = "Development"

import numpy as np
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
import sys, os, time

def main(args):
    try:
        fixFile = args[0]
        movFile = args[1]
        outfile = args[2]
        iterations = int(args[3]) if len(args) > 3 else 10000
        epsilon = float(args[4]) if len(args) > 4 else 0.1
    except Exception as e:
        print("Usage: ransac_2d_trafo.py fixFile.csv movFile.csv [iterations=10000] [epsilon=0.1]\n "
              "fixFile is a comma seperated file with X/Y coordinates of the fixed point cloud\n"
              "movFile is ------------------\"------------------------ of the moving point cloud")
        print("Error: %s" % e)

    fixPts = np.loadtxt(fixFile, delimiter=",")
    movPts = np.loadtxt(movFile, delimiter=",")
    num_movPts = movPts.shape[0]

    # create distance matrix for fixed points
    dists = distance_matrix(fixPts, fixPts)

    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axis('equal')
    fixplot, = ax.plot(fixPts[:, 0], fixPts[:, 1], 'bo')
    movplot, = ax.plot(movPts[:, 0], movPts[:, 1], 'rx')

    currplot, = ax.plot([], [], 'g*')
    fig.canvas.draw()

    curr_it = 0
    best_score = None
    best_solution = None
    while curr_it < iterations:
        curr_idx = np.random.choice(range(num_movPts), 2, replace=False)
        movp1 = movPts[curr_idx[0]]
        movp2 = movPts[curr_idx[1]]
        curr_dist = np.linalg.norm(movp2 - movp1)
        print("Curr_dist: %.3f" % curr_dist)
        curr_candidates = np.abs(dists - curr_dist) < epsilon
        print("Fitting pairs: %d" % np.count_nonzero(curr_candidates))
        xidxs, yidxs = np.nonzero(curr_candidates)
        # go over all fitting pairs
        for (xidx, yidx) in zip(xidxs, yidxs):
            fixp1 = fixPts[xidx]
            fixp2 = fixPts[yidx]
            # fit pair
            t = movp1 - fixp1
            vec_fix = fixp2 - fixp1
            vec_mov = movp2 - movp1
            rot = -np.arccos(np.inner(vec_mov, vec_fix) / (np.linalg.norm(vec_fix) * np.linalg.norm(vec_mov)))
            R = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])
            # transform all points
            movPts_tf = movPts - movp1
            movPts_tf = np.einsum('ij, jk -> ik', movPts_tf, R)
            movPts_tf = movPts_tf + fixp1
            # get score
            min_dists = np.min(distance_matrix(movPts_tf, fixPts), axis=0)
            #score = np.sum(np.square(min_dists[min_dists < epsilon]))
            score = np.sum(min_dists[min_dists < epsilon]) + np.count_nonzero(min_dists > epsilon)*1
            if best_score is None or score < best_score:
                print("Found a new best: score %.3f" % score)
                best_score = score
                best_solution = [t, R, movp1, fixp1]
                print('"Affine[%.8f %.8f]"' % (-movp1[0], -movp1[1]))
                print('"Affine[%.8f %.8f %.8f %.8f]"' % (R[0, 0], R[0,1], R[1,0], R[1,1]))
                print('"Affine[%.8f %.8f]"' % (fixp1[0], fixp1[1]))
                #print(best_solution)
                movplot.set_xdata(movPts_tf[:, 0])
                movplot.set_ydata(movPts_tf[:, 1])
                currplot.set_xdata([fixp1[0], fixp2[0]])
                currplot.set_ydata([fixp1[1], fixp2[1]])
                fig.canvas.draw()
                fig.canvas.flush_events()
                #time.sleep(0.500)
        curr_it += 1
    print("Done with %i iterations" % (iterations))
    plt.ioff()
    plt.show()

    outf = open(outfile, "w")
    outf.write('-filter "Affine[%.8f %.8f]"' % (-movp1[0], -movp1[1]))
    outf.write('-filter "Affine[%.8f %.8f %.8f %.8f]"' % (R[0, 0], R[0,1], R[1,0], R[1,1]))
    outf.wrtie('-filter "Affine[%.8f %.8f]"' % (fixp1[0], fixp1[1]))
    outf.close()

if __name__ == '__main__':
    main(sys.argv[1:])