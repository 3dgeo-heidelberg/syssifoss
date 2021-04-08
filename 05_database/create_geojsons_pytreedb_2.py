#!/usr/bin/env python
# -- coding: utf-8 --
# Hannah Weiser, Heidelberg University, February 2021
# h.weiser@stud.uni-heidelberg.de

"""
Script for generating geojsons for the SYSSIFOSS database from the folder, where individual trees and metrics are stored

Usage: create_geojsons_pytreedb.py <in_dir> <out_dir>
    in_dir      directory, in which individual tree files are stored
                - one folder per tree, named by the treeID in the format speciesID_plotID_number
                                e.g. AbiAlb_BR03_01
                - two position files per tree:
                    - WGS48: named as: speciesID_plotID_number_position.txt
                                e.g. AbiAlb_BR03_01_position.txt
                    - EPSG25832: named as speciesID_plotID_number_m_epsg25832.txt
                                AbiAlb_BR03_01_m_position_epsg25832.txt
                - multiple metric files, named as speciesID_plotID_number_m_YYYY-MM-DD_platform-season.txt,
                                e.g. AbiAlb_BR03_01_m_2019-08-24_ULS-on.txt
                                            or as speciesID_plotID_number_m_YYYY-MM-DD_FI.txt (for field measurements)
                                e.g. AbiAlb_BR03_01_m_2019-06-04_FI.txt
                - a "pointcloud" folder containing all available pointclouds,
                        named as speciesID_plotID_number_YYYY-MM-DD_quality_platform-season.laz
                                e.g. AbiAlb_BR03_01_2019-07-05_q3_ALS-on.laz
    out_dir     directory to which the GeoJSONS will be written (one file per tree)
"""

import glob
import os
import sys
import json
import laspy
import pandas as pd


def get_species_from_id(id):
    sp_id = id.split("_")[0]
    if sp_id == "FagSyl":
        return "Fagus sylvatica"
    elif sp_id == "CarBet":
        return "Carpinus betulus"
    elif sp_id == "BetPen":
        return "Betula pendula"
    elif sp_id == "PseMen":
        return "Pseudotsuga menziesii"
    elif sp_id == "QuePet":
        return "Quercus petraea"
    elif sp_id == "QueRob":
        return "Quercus robur"
    elif sp_id == "PicAbi":
        return "Picea abies"
    elif sp_id == "AceCam":
        return "Acer campestre"
    elif sp_id == "AcePse":
        return "Acer pseudoplatanus"
    elif sp_id == "PinSyl":
        return "Pinus sylvestris"
    elif sp_id == "LarDec":
        return "Larix decidua"
    elif sp_id == "QueRub":
        return "Quercus rubra"
    elif sp_id == "AbiAlb":
        return "Abies alba"
    elif sp_id == "FraExc":
        return "Fraxinus excelsior"
    elif sp_id == "PruAvi":
        return "Prunus avium"
    elif sp_id == "PruSer":
        return "Prunus serotina"
    elif sp_id == "JugReg":
        return "Juglans regia"
    elif sp_id == "SalCap":
        return "Salix caprea"
    elif sp_id == "TsuHet":
        return "Tsuga heterophylla"
    elif sp_id == "TilSpe":
        return "Tilia spec."
    elif sp_id == "SorTor":
        return "Sorbus torminalis"


def get_pcount_from_las(path):
    infile = laspy.file.File(path, mode="r")
    header = infile.header
    n_points = header.point_return_count

    return sum(n_points)


def get_pcount_from_xyz(path):
    df = pd.read_csv(path, header=0)

    return df.shape[0]


in_dir = sys.argv[1]
out_dir = sys.argv[2]

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

list_subfolders = [f.path for f in os.scandir(in_dir) if f.is_dir()]

for sub_dir in list_subfolders:
    os.chdir(sub_dir)
    full_tree_id = os.path.basename(sub_dir)
    id_parts = full_tree_id.split("_")
    outpath = os.path.join(out_dir, full_tree_id + ".geojson")

    if os.path.exists(outpath):
        continue

    if "broadl" in outpath:
        continue

    # basic structure of geojson
    tree_dict = {
        "type": "Feature",
        "properties": {
            "id": full_tree_id,
            "species": get_species_from_id(full_tree_id),
            "measurements": [

            ],
            "data": [

            ]
        },
        "geometry": {
            "type": "Point"
        }
    }
    season = None

    measurements = glob.glob("*_m*.txt")
    all_keys = []
    for i, entry in enumerate(measurements):
        sp_id, plot, tree_id, _, date, sensseas = os.path.basename(entry).split(".")[0].split("_")[:6]
        source = sensseas.split("-")[0]
        if len(sensseas.split("-")) > 1:
            season = "leaf-" + sensseas.split("-")[1]
        if "position" in entry:
            with open(os.path.join(sub_dir, entry), "r") as e:
                keys = ["position_xyz"]
                values = [[float(coord) for coord in e.readline().strip().split(",")]]
                tree_dict["properties"]["measurements"].append({
                    "crs": "epsg:25832"
                })
        else:
            with open(os.path.join(sub_dir, entry), "r") as e:
                keys = e.readline().strip().split()
                values = e.readline().strip().split()
                if season is not None:
                    tree_dict["properties"]["measurements"].append({
                        "source": source,
                        "date": date,
                        "canopy_condition": season
                    })
                elif season is None:
                    tree_dict["properties"]["measurements"].append({
                        "source": source,
                        "date": date,
                    })
        for j, key in enumerate(keys):
            tree_dict["properties"]["measurements"][i][key] = values[j]
            all_keys.append(key)
    if "height_m" not in all_keys:
        print("Required metrics missing. Skipping %s." % full_tree_id)
        continue

    if glob.glob("*position.txt"):
        pos_file = glob.glob("*position.txt")[0]
    else:
        print("Position file missing. Skipping %s." % full_tree_id)
        continue
    with open(os.path.join(sub_dir, pos_file), "r") as p:
        x, y, z = p.readline().strip().split(",")
        tree_dict["geometry"]["coordinates"] = [float(coord) for coord in [y, x, z]]

    data_folders = [f.path for f in os.scandir(sub_dir) if f.is_dir()]
    for data_folder in data_folders:
        datatype = data_folder.split("\\")[-1]
        entries = [f.name for f in os.scandir(data_folder) if f.is_file()]
        if len(entries) > 0 and datatype == "pointcloud":
            for entry in entries:
                if "READ_ME" in entry:
                    continue
                sp_id, plot, tree_id, date, qual, sensseas = os.path.basename(entry).split(".")[0].split("_")
                mode = sensseas.split("-")[0]
                if mode == "ALS":
                    sensor = "RIEGL LMS-VQ780i"
                elif mode == "ULS":
                    sensor = "RIEGL miniVUX-1UAV"
                elif mode == "TLS":
                    sensor = "RIEGL VZ-400"
                quality = qual[1]
                if os.path.splitext(entry)[1] == ".laz":
                    p_count = get_pcount_from_las(os.path.join(data_folder, entry))
                season = "leaf-" + sensseas.split("-")[1]
                tree_dict['properties']["data"].append({
                    "type": datatype,
                    "mode": mode,
                    "date": date,
                    "canopy_condition": season,
                    "sensor": sensor,
                    "point_count": p_count,
                    "crs": "epsg:25832",
                    "file": entry,
                    "quality": quality
                })

    print("Writing %s..." % outpath)
    with open(outpath, 'w') as outfile:
        output = json.dumps(tree_dict, indent=4)
        outfile.write(output)
