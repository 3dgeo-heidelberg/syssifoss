#!/usr/bin/env python
# -- coding: utf-8 --
# Hannah Weiser, Heidelberg University, February 2021
# h.weiser@stud.uni-heidelberg.de

"""
Script to write tree properties from SYSSIFOSS database GeoJSONs to tab-delimited csv-files.
One file is generated for general properties (species and position in different coordinate systems)
                    + for specific tree metrics from different sources (i.e. point clouds or field measurements)

Usage:
geojsons_to_csv.py inf
inf         input geojson files (wildcards can be used to specify multiple files)

Output CSV-files are written to the location of the input files.
"""

import glob
import sys
import json
import pandas as pd
import numpy as np

inf = sys.argv[1]
infiles = glob.glob(inf)

for inf in infiles:
    # read file
    with open(inf, 'r') as infile:
        data = infile.read()

    # parse file
    tree_dict = json.loads(data)

    # outfiles
    outf_general = inf.replace(".geojson", "_general.txt")
    outf_metrics = inf.replace(".geojson", "_metrics.txt")

    df_general = pd.DataFrame()
    df_general["species"] = [tree_dict["properties"]["species"]]
    df_general["lat_epsg4326"] = [tree_dict["geometry"]["coordinates"][1]]
    df_general["long_epsg4326"] = [tree_dict["geometry"]["coordinates"][0]]
    df_general["elev_epsg4326"] = [tree_dict["geometry"]["coordinates"][2]]

    n = len(tree_dict["properties"]["measurements"])
    df_metrics = pd.DataFrame(index=np.arange(n-1))
    i = 0
    for entry in tree_dict["properties"]["measurements"]:
        columns = entry.keys()
        if "position_xyz" in columns:
            df_general["x_epsg25832"] = entry["position_xyz"][0]
            df_general["y_epsg25832"] = entry["position_xyz"][1]
            df_general["z_epsg25832"] = entry["position_xyz"][2]
        else:
            for col in columns:
                df_metrics.loc[i, col] = entry[col]
            i += 1

    # write outfiles
    df_general.to_csv(outf_general, index=False, sep="\t")
    df_metrics.replace("NA", "", inplace=True)
    df_metrics.to_csv(outf_metrics, index=False, sep="\t")

