# SYSSIFOSS

SYSSIFOSS is a joint project between the [Institute of Geography and Geoecology (IFGG) of the Karlsruhe Institute of Technology (KIT)](https://www.ifgg.kit.edu/english/index.php) and the [3DGeo Research Group of Heidelberg University](https://www.geog.uni-heidelberg.de/3dgeo/index_en.html). In this project we suggest a new approach to create synthetic LiDAR data by combining the outputs of an established forest growth simulator with a database of species-specific model trees extracted from real LiDAR point clouds.

Find out more on our [project website](www.uni-heidelberg.de/syssifoss) and watch our [project video](https://www.youtube.com/watch?v=B9yStyUBaa0&t=23s).

The repository contains the scripts used for processing our laser scanning forest data. 

## Repository structure:

### 01_tls_georeferencing
- **OPALS** tools for registering TLS data to ULS data: 
    - stem detection: `opals_stem_detection_tls.py`, `opals_stem_detection_uls.py`      
    - applying multiple transformation: `opals_apply_transformation.py`
      
    - running the ICP: `opals_run_icp.py` and `opals_run_icp_all.py`
- Execute these scripts in the **OPALS Shell** or from the PyCharm-distribution that comes with Opals
- RANSAC script for 2D stem matching: `ransac_2d_trafo.py`

### 02_tree_extraction:
- Scripts for automatically extracting point clouds from one point cloud using template point clouds: 
  `get2clouds_multiple_targets.py`, `get2clouds_ULS_ALS.py`
Can be used for
    - retrieving all points of a spatial subset from a large point cloud using a downsampled source (template) point cloud from the same dataset,
    e.g.,
        - *template*: Tree point cloud segmented from downsampled TLS point cloud (only x, y, z)
        - *source*: Full TLS point cloud (x, y, z, + other attributes like reflectance, return number, etc.)
        - *output*: Tull-resolution TLS tree point cloud (x, y, z + other attributes)
    
    - retrieving all points of a spatial subset from one dataset using a source point cloud from another dataset, e.g.,
        - *template*: Tree point cloud segmented from TLS point cloud
        - *source*: ULS point cloud of a full forest plot
        - *output*: ULS single tree point cloud
    
### 03_metric_computation:
- Library for the computation of tree metrics: `TreeMetrics.py`

- Scripts for 
    - computing tree metrics, using `TreeMetrics.py`: `compute_tree_metrics.py`
      
    - computating of DBH: `DBH_RANSAC.py`, `DBH_ellipse.py`
    
    - deriving tree positions: `tree_positiony.py`
    
- External programs used:
    - Concave hull C++ program by Alasdair Craig: https://www.codeproject.com/Articles/1201438/The-Concave-Hull-of-a-Set-of-Points
    using the algorithm by Moreira & Santos (2007): `concave.exe`
      
    - Ellipse fitting program by Nicky van Foreest:https://github.com/ndvanforeest/fit_ellipse: `fit_ellipse.py`

### 04_opals_formatdef:
- OPALS format definitions specific to the SYSSIFOSS data.


### 05_database
- Scripts to write GeoJSON files and CSV files for the database from the contents of a tree data folder: 
    - GeoJSONs: `create_geojsons_pytreedb.py` and `create_geojsons_pytreedb_2.py`
    - CSV: `geojsons_to_csv.py`

## Requirements

Requirements for all non-opals scripts are found in [requirements.txt](requirements.txt)

## License

See [LICENSE](LICENSE)