# Prepare packages -----------------------------------------------------------------------------------------------------
# Packages for paths manipulation
from folder import *                              # For paths of sub-folders
import os                                         # For manipulating the directory path and to execute commands
import pathlib                                    # For manipulating the directory path

# Packages for untransformation
import numpy as np                                # For all calculation and data array/matrices manipulation
from transformation import center_calculation     # For center calculation

# Packages for unrotating and untranslating
import rasterio.features                          # For vectorising features in array
from shapely.geometry import shape                # For manipulating spatial information (geometry) under GeoJSON format
import geopandas as gpd                           # For manipulating shape files
from pyogrio import write_dataframe               # For writing out shape file (twice faster than geopandas)

# Package for manipulating spatial data
import rasterio                                   # For reading and manipulating spatial data

# ----------------------------------------------------------------------------------------------------------------------


def polygons_untransformation(transformation_selection, number_simulation, angle_func, x_translation_func,
                              y_translation_func, time_extract_func):
    """This function is to unrotate and write rasters into shapefiles

    -----------
    References: https://docs.dea.ga.gov.au/notebooks/Frequently_used_code/Polygonise_pixel_edges.html
                https://gis.stackexchange.com/questions/187877/how-to-polygonize-raster-to-shapely-polygons/

                https://github.com/sgillies/affine/blob/master/affine/__init__.py#L178
                https://gis.stackexchange.com/questions/408386/rotating-raster-using-python
                https://pythonhosted.org/PyAgg/affine.m.html

                https://gis.stackexchange.com/questions/408386/rotating-raster-using-python
                https://gis.stackexchange.com/questions/350526/gdal-setgeotransform-issue
                https://corteva.github.io/rioxarray/stable/examples/convert_to_raster.html
                https://gdal.org/tutorials/geotransforms_tut.html#geotransforms-tut
                https://gdal.org/tutorials/raster_api_tut.html
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                number_simulation:
                (int)
                                                Ordinal number of simulation
                angle_func:
                (int or float)
                                                Angle that is used to transform LiDAR data
                x_translation_func:
                (int or float)
                                                x coordinate that is used to transform LiDAR data
                y_translation_func:
                (int or float)
                                                y coordinate that is used to transform LiDAR data
                time_extract_func:
                (int)
                                                Amount of time that flood model predicted
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed = "rotated"
        untransformed_path = unrotated_path
        flowdepth_path = rotated_flowdepth
    elif transformation_selection == 't':
        transformed = "translated"
        untransformed_path = untranslated_path
        flowdepth_path = translated_flowdepth
    else:
        transformed = "combined"
        untransformed_path = uncombined_path
        flowdepth_path = combined_flowdepth

    # Extract the array from COMBINED raster
    raster_poly = rasterio.open(fr"{flowdepth_path}\\flowdepth_{transformed}_{number_simulation}_at_{time_extract_func}.nc")
    raster_array = raster_poly.read(1)
    raster_transform = raster_poly.transform
    raster_crs = raster_poly.crs

    # Extract parameters: id and depth
    id_pixels = np.arange(raster_array.size).reshape(raster_array.shape)
    depth_list = list(raster_array.flatten())

    # UNCOMBINE 1 (untranslate) affine transformation
    uncombined1_transform = raster_transform * raster_transform.translation(x_translation_func, y_translation_func)

    # Find out centers
    raster_center = center_calculation(transformation_selection, False)

    # UNCOMBINE 2 (unrotate) affine transformation
    uncombined2_transform = uncombined1_transform * uncombined1_transform.rotation(angle_func, raster_center)

    # Change the name
    uncombined_transform = uncombined2_transform

    # Vectorise features
    uncombined_vectors = rasterio.features.shapes(source=id_pixels.astype(np.int16),
                                                  transform=uncombined_transform)

    # List the UNCOMBINED polygons to extract necessary parameters
    uncombined_vectors_list = list(uncombined_vectors)

    # Get geometry
    polygons_geometry_values = [shape(polygon_geometry) for polygon_geometry, value_geometry in uncombined_vectors_list]

    # Get id
    id_uncombined_pixels_values = [id_value for id_polygon, id_value in uncombined_vectors_list]

    # Create UNCOMBINED database under geopandas dataframe (gdf) format
    uncombined_data = {"id": id_uncombined_pixels_values,
                       "depth": depth_list}
    uncombined_raster_poly_gdf = gpd.GeoDataFrame(data=uncombined_data,
                                                  geometry=polygons_geometry_values,
                                                  crs=raster_crs)

    shapefile_dir = pathlib.Path(os.getcwd()) / pathlib.Path(f"{untransformed_path}\\shapefile_{number_simulation}")
    if not os.path.exists(shapefile_dir):
        os.mkdir(shapefile_dir)

    write_dataframe(uncombined_raster_poly_gdf,
                    fr"{untransformed_path}\\shapefile_{number_simulation}\\flowdepth_un{transformed}_{number_simulation}_at_{time_extract_func}.geojson",
                    driver="GeoJSON")