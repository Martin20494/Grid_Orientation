# Prepare packages -----------------------------------------------------------------------------------------------------
# For paths manipulation
from folder import *                                      # For paths of sub-folders
import os                                                 # For creating new folders
import pathlib                                            # For getting path

# For handling polygons
import geopandas as gpd                                   # For manipulating geospatial data under polygon format
from shapely.geometry import Point                        # For converting point coordinates into shapely geometry point

# For clipping
from datasetClipping import clip, xyz_array               # For clipping dataset

# ----------------------------------------------------------------------------------------------------------------------

# Create a function to convert raster array into shape file
def raster_array_to_shapefile(transformation_selection,
                              number_simulation,
                              time_extract_func,
                              adjusted_value_list=[0, 0, 0, 0],
                              clip_using_shape_file=False):
    """This function is to convert a raster array into a shape file

    -----------
    References: https://geopandas.org/docs/reference/api/geopandas.GeoDataFrame.to_crs.html
                https://shapely.readthedocs.io/en/stable/manual.html
                https://geopandas.org/docs/user_guide/io.html
                https://sgillies.net/2014/01/18/getting-shapes-of-raster-features-with-rasterio.html
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
                time_extract_func:
                (int)
                                            Amount of time that flood model predicted
                adjusted_value_list:
                (list)
                                            List contain values to change the boundaries
                                            [distance xmin, distance ymin, distance xmax, distance ymax]
                clip_using_shape_file:
                (boolean)
                                            Selection for clipping by using shapfile or manually
                                            True means using shapefile and clipping out all -999/nodata values
                                            False means clipping manually (following rectangle shape)
                                            Default is False
    -----------

    -----------
    Return:
               clipped_raster_func:
                                            A dataset without padding
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        untransformed_path = unrotated_path
    elif transformation_selection == 't':
        untransformed_path = untranslated_path
    else:
        untransformed_path = uncombined_path

    # Call out and clip the raster array
    raster_func = xyz_array(transformation_selection, number_simulation, time_extract_func)
    clipped_raster_func = clip(transformation_selection, raster_func, adjusted_value_list, clip_using_shape_file)

    # Convert x, y coordinates array into shapely geometry
    point_geo_values_func = [Point(clipped_raster_func[i, 0], clipped_raster_func[i, 1]) for i in
                             range(clipped_raster_func.shape[0])]

    # Build up geopandas dataframe
    point_data_func = {"depth": clipped_raster_func[:, 2]}
    point_gdf = gpd.GeoDataFrame(data=point_data_func,
                                 geometry=point_geo_values_func,
                                 crs=2193)

    shapefile_dir = pathlib.Path(os.getcwd()) / pathlib.Path(f"{untransformed_path}\\shapefile_Point")
    if not os.path.exists(shapefile_dir):
        os.mkdir(shapefile_dir)

    # Write geopandas dataframe into shape file
    point_gdf.to_file(f"{untransformed_path}\\shapefile_Point\\Point_clip.shp",
                      driver="ESRI Shapefile")

    return clipped_raster_func