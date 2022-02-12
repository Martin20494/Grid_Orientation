# Prepare packages -----------------------------------------------------------------------------------------------------
# Main packages for downloading lidar
import geoapis
import geoapis.lidar

# Other packages
from folder import *  # For paths of sub-folders
import os  # For manipulating the directory path
import pathlib  # For creating and manipulating the directory path
import geopandas  # For creating a series to store geometry object
import shapely.geometry  # For creating a polygon object
import shutil  # For copying file/folder
# ----------------------------------------------------------------------------------------------------------------------


def download_lidar(boundary_coordinates_func, name_dataset, lidar_number):
    """This function is to download las/laz files from Open Topography

    -----------
    References:  https://github.com/niwa/geoapis/wiki/Basic-Usage
    -----------

    -----------
    Arguments:
                boundary_coordinates_func:
                (list)
                                            A list of xmin, ymin, xmax, ymax of rectangular boundary
                name_dataset:
                (string)
                                            Name of the dataset
                lidar_number:
                (int)
                                            Ordinal number of lidar folder
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Get crs
    h_crs = 2193

    # Get xmin, ymin, xmax, ymax of rectangular boundary
    x0_func = boundary_coordinates_func[0]
    y0_func = boundary_coordinates_func[1]
    x1_func = boundary_coordinates_func[2]
    y1_func = boundary_coordinates_func[3]

    # Assign those coordinates into shapely geometry polygon and store under geopandas format
    lidar_bound_coordinates = shapely.geometry.Polygon(
        [(x0_func, y0_func), (x1_func, y0_func), (x1_func, y1_func), (x0_func, y1_func)])
    lidar_bound_coordinates = geopandas.GeoSeries([lidar_bound_coordinates])
    lidar_bound_coordinates = lidar_bound_coordinates.set_crs(h_crs)

    # Create lidar folder to store necessary files
    lidar_folder_dir = pathlib.Path(os.getcwd()) / pathlib.Path(f"{original_lidar_path}\\lidar_{lidar_number}")
    if not os.path.exists(lidar_folder_dir):
        os.mkdir(lidar_folder_dir)

    # Create test file to store test_catchment.zip file (one of necessary files for downloading las/laz files)
    test_dir = pathlib.Path(os.getcwd()) / pathlib.Path(f"{original_lidar_path}\\lidar_{lidar_number}\\test")
    if not os.path.exists(test_dir):
        os.mkdir(test_dir)

    # Create test_catchment.zip
    test_path = pathlib.Path(f"{original_lidar_path}\\lidar_{lidar_number}\\test\\test_catchment")
    lidar_bound_coordinates.to_file(test_path)
    shutil.make_archive(base_name=test_path, format='zip', root_dir=test_path)
    shutil.rmtree(test_path)
    lidar_bound_path = pathlib.Path(str(test_path) + ".zip")

    # Create polygon
    lidar_polygon = geopandas.read_file(lidar_bound_path)
    lidar_polygon.to_crs(h_crs)

    # Download lidar
    lidar_fetcher = geoapis.lidar.OpenTopography(cache_path=fr"{original_lidar_path}\\lidar_{lidar_number}",
                                                 search_polygon=lidar_polygon,
                                                 verbose=True)
    lidar_fetcher.run(name_dataset)
