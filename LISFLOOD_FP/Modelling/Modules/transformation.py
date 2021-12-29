# Prepare packages -----------------------------------------------------------------------------------------------------
# Package for paths manipulation
from folder import *                                # For paths of sub-folders
import glob                                         # For getting all files' names in a certain path
import os                                           # For manipulating paths and to execute commands
import pathlib                                      # For manipulating the directory path
import shutil                                       # For copying file/folder

# Package for reading las/laz file and handling lidar point cloud
import numpy as np                                  # For all calculation and data array/matrices manipulation
import pdal                                         # For handling LiDAR point cloud data
import json                                         # For creating json text data format

# Packages for transformation
from numba import guvectorize, float64              # For speeding up the time of the code running

# Packages for unrotating and untranslating
from osgeo import gdal                              # For manipulating rasters (calculating centers)
from shapely.geometry import Polygon                # For creating polygons to change raster values
import geopandas as gpd                             # For manipulating shape files

# Package for timing steps
import time                                         # For timing steps
# ----------------------------------------------------------------------------------------------------------------------


# Remove all warnings --------------------------------------------------------------------------------------------------
# Reference: https://numba.pydata.org/numba-doc/dev/reference/deprecation.html#suppressing-deprecation-warnings
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaWarning)
# ----------------------------------------------------------------------------------------------------------------------



# READING LAS/LAZ FILES OF TILES
def read_las_file(filename_tile):
    """ This function is to read the las/las file into arrays

    -----------
    References: https://pythonhosted.org/laspy/file.html
                https://pythonhosted.org/laspy/tut_part_1.html
    -----------

    -----------
    Arguments:
                filename_tile:
                (string)
                                                Name of las/laz files of tiles

    Returns:
                points_classification_func:
                (array)
                                                An array includes classification information (int format)
                coordinates_func:
                (array)
                                                An array includes x, y, z coordinates values of points as rows
    """
    # Set up horizontal and vertical crs
    h_crs = 2193
    v_crs = 7839

    # Write las/las file into arrays
    pdal_pipeline_instructions = [{"type": "readers.las", "filename": fr"{filename_tile}"},
                                  {"type": "filters.reprojection", "in_srs": f"EPSG:{h_crs}+{v_crs}",
                                   "out_srs": f"EPSG:{h_crs}+{v_crs}"}]

    pdal_pipeline = pdal.Pipeline(json.dumps(pdal_pipeline_instructions))
    pdal_pipeline.execute()

    # Call out arrays
    points_full_func = pdal_pipeline.arrays  # note is shape [1,nx,ny]

    # Assign x, y , z (values in meters)
    x_func = points_full_func[0]["X"]
    y_func = points_full_func[0]["Y"]
    z_func = points_full_func[0]["Z"]

    # Get classification
    point_classification_func = points_full_func[0]['Classification']

    # Write x, y, z into arrays with each row as (x, y, z) (meters as unit)
    coordinates_func = np.vstack((x_func, y_func, z_func)).transpose().astype('float64')

    return point_classification_func, coordinates_func


# CENTER CALCULATION
def center_calculation(transformation_selection, lidar=True):
    """This function is to calculate the center of reference DEM raster which will be used in the whole process

    -----------
    Reference:  https://gis.stackexchange.com/questions/104362/how-to-get-extent-out-of-geotiff
                https://rasterio.readthedocs.io/en/latest/quickstart.html
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                lidar:
                (boolean)
                                            True means doing rotation on LiDAR data
                                            False means doing unrotation on model outputs
    -----------

    -----------
    Returns:
                (tuple)
                                            A tuple of x, y coordinates of center raster
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed_nc_raster_path = rotated_nc_raster_path
    elif transformation_selection == 't':
        transformed_nc_raster_path = translated_nc_raster_path
    else:
        transformed_nc_raster_path = combined_nc_raster_path

    # Read the file using gdal
    raster_reference_func = gdal.Open(fr"{transformed_nc_raster_path}\\generated_dem_reference.nc")

    if lidar:
        # Get the center coordinates for LiDAR
        xmin, xpixel, _, ymax, _, ypixel = raster_reference_func.GetGeoTransform()
        width, height = raster_reference_func.RasterXSize, raster_reference_func.RasterYSize
        xmax = xmin + width * xpixel
        ymin = ymax + height * ypixel
        center_x_func = ((xmin + xmax) / 2)
        center_y_func = ((ymin + ymax) / 2)
        center_func = np.array([center_x_func, center_y_func])

    else:
        # Get the center coordinates for DEMs
        center_func = (raster_reference_func.RasterXSize / 2, raster_reference_func.RasterYSize / 2)
    return center_func


# ROTATION #############################################################################################################
# Create a command for guvectorize() decorator
# In that, float64[:,:] stands for the format of the array and float64 stands for the format of other parameters
# (m,n) signature stands for the matrix of the array and () signature stands for other parameters
gu_rotation = guvectorize([(float64[:, :], float64, float64, float64, float64, float64[:, :])],
                          '(m,n),(),(),(),()->(m,n)')

def point_rotation(coordinates_func, angle, center_x_func, center_y_func, clockwise, new_coordinates_func):
    """This function is to calculate the rotated coordinates of lidar data

    -----------
    References: https://www.youtube.com/watch?v=RqZH-7hlI48
                https://stackoverflow.com/questions/14607640/rotating-a-vector-in-3d-space
                https://stackoverflow.com/questions/5954603/transposing-a-1d-numpy-array
                https://en.wikipedia.org/wiki/Rotation_matrix
                https://math.stackexchange.com/questions/270194/how-to-find-the-vertices-angle-after-rotation

                https://github.com/numba/numba/issues/3312
                https://numba.pydata.org/numba-doc/latest/user/vectorize.html
                https://numba.pydata.org/numba-doc/latest/cuda/ufunc.html
                http://numba.pydata.org/numba-doc/0.20.0/reference/compilation.html
                http://numba.pydata.org/numba-doc/0.12/tutorial_numpy_and_numba.html
                https://numba.pydata.org/numba-doc/dev/reference/types.html
    -----------

    -----------
    Arguments:
                coordinates_func:
                (array)
                                        An array of the coordinates of the point will be rotated
                angle:
                (int or float)
                                        The value of angle to rotate
                center_x_func:
                (float)
                                        Coordinate value of x center.
                                        The x center here was used from the center of reference DEM without padding
                center_y_func:
                (float)
                                        Coordinate value of x center.
                                        The x center here was used from the center of reference DEM without padding
                clockwise:
                (string)
                                        Rotating the points in clockwise (1) or anti_clockwise (0) directions
                new_coordinates_func:
                (array)
                                        A new array of rotated x, y, z coordinates values of a point
    -----------

    -----------
    Returns:
                (array)
                                        The function will return itself without the need of 'return' command
    -----------

    """
    # Convert degree to radian and calculate cosine and sine
    radian = np.deg2rad(angle)
    cosine = np.cos(radian)
    sine = np.sin(radian)

    # Create a for loop to manipulate each row of the array
    for i in range(coordinates_func.shape[0]):
        # Do subtraction with center coordinates
        diff_x = coordinates_func[i, 0] - center_x_func
        diff_y = coordinates_func[i, 1] - center_y_func

        # Calculate the rotated point coordinates
        if clockwise == 0:  # Rotating in anti-clockwise direction
            new_coordinates_func[i, 0] = diff_x * cosine - diff_y * sine + center_x_func
            new_coordinates_func[i, 1] = diff_x * sine + diff_y * cosine + center_y_func
            new_coordinates_func[i, 2] = coordinates_func[i, 2]
        else:  # Rotating in clockwise direction
            new_coordinates_func[i, 0] = diff_x * cosine + diff_y * sine + center_x_func
            new_coordinates_func[i, 1] = diff_x * (-sine) + diff_y * cosine + center_y_func
            new_coordinates_func[i, 2] = coordinates_func[i, 2]


# Wrapping function to map later
wrapping_point_rotation = gu_rotation(point_rotation)

# END ROTATION #########################################################################################################


# TRANSLATION ##########################################################################################################
# Create a command for guvectorize() decorator
# In that, float64[:,:] stands for the format of the array and float64 stands for the format of other parameters
# (m,n) signature stands for the matrix of the array and () signature stands for other parameters
gu_translation = guvectorize([(float64[:, :], float64, float64, float64[:, :])], '(m,n),(),()->(m,n)')

# Create a point translation function
def point_translation(coordinates_func, x_translation_func, y_translation_func, new_coordinates_func):
    """This function is to calculate the translated coordinates of lidar data

    -----------
    References: https://www.youtube.com/watch?v=RqZH-7hlI48
                https://stackoverflow.com/questions/14607640/rotating-a-vector-in-3d-space
                https://stackoverflow.com/questions/5954603/transposing-a-1d-numpy-array
                https://en.wikipedia.org/wiki/Rotation_matrix
                https://math.stackexchange.com/questions/270194/how-to-find-the-vertices-angle-after-rotation

                https://github.com/numba/numba/issues/3312
                https://numba.pydata.org/numba-doc/latest/user/vectorize.html
                https://numba.pydata.org/numba-doc/latest/cuda/ufunc.html
                http://numba.pydata.org/numba-doc/0.20.0/reference/compilation.html
                http://numba.pydata.org/numba-doc/0.12/tutorial_numpy_and_numba.html
                https://numba.pydata.org/numba-doc/dev/reference/types.html
    -----------

    -----------
    Arguments:
                coordinates_func:
                (array)
                                        An array of the coordinates of the point will be rotated
                x_translation_func:
                (int or float)
                                        Distances to translate x coordinates values
                y_translation_func:
                (int or float)
                                        Distances to translate y coordinates values
                new_coordinates_func:
                (array)
                                        A new array of translated x, y, z coordinates values of a point
    -----------

    -----------
    Returns:
               (array)
                                        The function will return itself without the need of 'return' command
    -----------

    """
    for i in range(coordinates_func.shape[0]):
        new_coordinates_func[i, 0] = coordinates_func[i, 0] + x_translation_func
        new_coordinates_func[i, 1] = coordinates_func[i, 1] + y_translation_func
        new_coordinates_func[i, 2] = coordinates_func[i, 2]

# Wrapping function to map later
wrapping_point_translation = gu_translation(point_translation)

# END TRANSLATION ######################################################################################################


# WRITING TRANSFORMED LAS/LAZ FILES ####################################################################################
def writing_las_file(transformation_selection, point_classification_func, number_simulation, transformed_array,
                     filename_code):
    """This function is to create a new lidar las/laz file
    the same as the old lidar las/laz file with different coordinates

    -----------
    References: https://pdal.io/tutorial/las.html
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                point_classification_func:
                (array)
                                            An array contains points' information from original tiles
                number_simulation:
                (int)
                                            Ordinal number of simulation
                transformed_array:
                (array)
                                            A transformed array
                filename_code:
                (string)
                                            Code of file name (the same as original lidar tiles' names)
    -----------

    -----------
    Returns:
                None.
    -----------
    """
    # Set up horizontal and vertical crs
    h_crs = 2193
    v_crs = 7839

    # Output file
    if transformation_selection == 'r':
        pathlib.Path(fr"{rotated_lidar_path}\\rotated_lidar_{number_simulation}\\Wellington_2013")\
            .mkdir(parents=True, exist_ok=True)
        outfile = fr"{rotated_lidar_path}\\rotated_lidar_{number_simulation}\\Wellington_2013\\{filename_code}.laz"
    elif transformation_selection == 't':
        pathlib.Path(fr"{translated_lidar_path}\\translated_lidar_{number_simulation}\\Wellington_2013")\
            .mkdir(parents=True, exist_ok=True)
        outfile = fr"{translated_lidar_path}\\translated_lidar_{number_simulation}\\Wellington_2013\\{filename_code}.laz"
    else:
        pathlib.Path(fr"{combined_lidar_path}\\combined_lidar_{number_simulation}\\Wellington_2013")\
            .mkdir(parents=True, exist_ok=True)
        outfile = fr"{combined_lidar_path}\\combined_lidar_{number_simulation}\\Wellington_2013\\{filename_code}.laz"

    # Get the shape of transformed LiDAR array
    lidar_arr_shape = len(point_classification_func)

    # Create empty transformed LiDAR array
    lidar_arr = np.empty([lidar_arr_shape], dtype=[('X', '<f8'), ('Y', '<f8'), ('Z', '<f8'), ('Classification', 'u1')])

    # Overwrite the X, Y, and Z values in the points array with the transformed values
    lidar_arr['X'] = (transformed_array[:, 0]).flatten()
    lidar_arr['Y'] = (transformed_array[:, 1]).flatten()
    lidar_arr['Z'] = (transformed_array[:, 2]).flatten()
    lidar_arr['Classification'] = point_classification_func

    # Save the file
    pdal_pipeline_instructions = [
        {
            "type": "writers.las",
            "scale_x": 0.0001,
            "scale_y": 0.0001,
            "offset_x": 1700000,
            "offset_y": 5400000,
            "a_srs": f"EPSG:{h_crs}+{v_crs}",
            "filename": str(outfile),
            "compression": "laszip"
        }
    ]

    # Run
    pdal_pipeline = pdal.Pipeline(json.dumps(pdal_pipeline_instructions), [lidar_arr])
    pdal_pipeline.execute()

# END WRITING ##########################################################################################################


# OTHER NECESSARY FUNCTIONS ############################################################################################
def tile_index_polygon(tile_array_coordinate):
    """This function is to create transformed tile index polygon for each original tile

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                tile_array_coordinate:
                (array)
                                        An array of each tile's coordinate (x, y, z)
    -----------

    -----------
    Returns:
                 tile_boundary_polygon:
                 (polygon)
                                        A polygon represents the boundary of the LiDAR tiles
    -----------

    """
    # Get xmin, xmax, ymin, ymax for the index
    x_min_tile = np.min(tile_array_coordinate[:, 0])
    x_max_tile = np.max(tile_array_coordinate[:, 0])
    y_min_tile = np.min(tile_array_coordinate[:, 1])
    y_max_tile = np.max(tile_array_coordinate[:, 1])

    # Design boundary coordinates
    tile_boundary = np.array([
        [x_min_tile, y_max_tile],
        [x_max_tile, y_max_tile],
        [x_max_tile, y_min_tile],
        [x_min_tile, y_min_tile],
        [x_min_tile, y_max_tile]
    ]).astype('float64')

    # Convert to polygon
    tile_boundary_polygon = Polygon(tile_boundary)

    return tile_boundary_polygon


def get_url(lidar_number):
    """This function is to extract the URL of each tile from original tile index file

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                lidar_number:
                (int)
                                Ordinal number of lidar folder
    -----------

    -----------
    Returns:
                url_file:
                (list)
                                URL of the specific tile
    -----------

    """
    # Read the shape file to have geopandas dataframe
    full_file = gpd.read_file(
        fr"{original_lidar_path}\\lidar_{lidar_number}\\Wellington_2013\\Wellington_2013_TileIndex.zip")

    # Get list of tiles names
    list_tile_paths = glob.glob(fr"{original_lidar_path}\\lidar_{lidar_number}\\Wellington_2013\*.laz")
    length_list = len(list_tile_paths)
    list_tile_names = [pathlib.Path(list_tile_paths[name]).stem for name in range(length_list)]

    url_file = []

    # Collect point classification and coordinate array
    for each_tile in range(length_list):
        # Get the row with the file name the same as given 'filename_tile'
        file_info = full_file[full_file['Filename'] == f"{list_tile_names[each_tile]}.laz"]

        # Get the index of that row
        url_index = file_info.index[0]

        # Use the extracted index to get the correct URL
        url_file.append(full_file['URL'][url_index])

    return url_file



def get_dictionary(lidar_number, type_dictionary):
    """This function is to create transformed tiles

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                lidar_number:
                (int)
                                                Ordinal number of lidar folder
                type_dictionary:
                (string)
                                                There are two types - classification and coordinate
    -----------

    -----------
    Return:
                point_dict:
                (dictionary)
                                                A dictionary of tiles about classification or coordinate
    -----------

    """
    # Get list of tiles names
    list_tile_paths = glob.glob(fr"{original_lidar_path}\\lidar_{lidar_number}\\Wellington_2013\*.laz")
    length_list = len(list_tile_paths)
    list_tile_names = [pathlib.Path(list_tile_paths[name]).stem for name in range(length_list)]

    # Get dictionary
    dict_point = dict()

    # Print big titles
    if type_dictionary == "classification":
        print("*** CLASSIFICATION DICTIONARY", "-" * 50)
    else:
        print("*** COORDINATE DICTIONARY", "-" * 50)

    # Start timing the whole processing
    whole_start = time.time()

    for each_tile in range(length_list):
        # Start timing of reading
        read_start = time.time()
        # Read each tile
        points_tiles, coordinates_tiles = read_las_file(list_tile_paths[each_tile])
        # End timing of reading
        read_end = time.time()

        # Select classification or coordinate
        if type_dictionary == "classification":
            # Add values in classification dictionary
            dict_point[f'{list_tile_names[each_tile]}.laz'] = points_tiles

        else:
            # Add values in coordinate dictionary
            dict_point[f'{list_tile_names[each_tile]}.laz'] = coordinates_tiles

        # PRINT -----------------------------------------------------------------------------------------
        print("*", "Tile", each_tile, "-", list_tile_names[each_tile], "-" * 30)

        # Get the running time of creating dictionaries
        read_time = read_end - read_start

        # Convert time to minute
        read_minute = read_time / 60

        # Time title
        read_title = "* Running time of creating dictionaries:"

        # Print
        print('{0:<60}{1:>10}'.format(read_title, read_minute))
        print("-" * 30)

    # End timing the whole processing
    whole_end = time.time()

    # Print --------------------------------------------------------------------------------------------
    # Get time of whole processing
    whole_time = whole_end - whole_start

    # Convert to minute
    whole_time_minute = whole_time / 60

    # Get name
    whole_name = "* Running time of whole dictionary generation process:"

    # Print
    print('{0:<60}{1:>10}'.format(whole_name, whole_time_minute))
    print("=" * 90)
    print("\n")

    return dict_point

# END OTHER NECESSARY FUNCTIONS ########################################################################################



# TRANSFORMATION EXECUTION #############################################################################################
def transformed_tiles_generation(transformation_selection, lidar_number,
                                 angle_func, x_translation_func, y_translation_func,
                                 center_x_func, center_y_func,
                                 tile_url_list,
                                 point_classification_func, point_coordinate_func,
                                 number_simulation):
    """This function is to create transformed tiles

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                lidar_number:
                (int)
                                                Ordinal number of lidar folder
                angle_func:
                (int or float)
                                                Angle that is used to transform LiDAR data
                x_translation_func:
                (int or float)
                                                x coordinate that is used to transform LiDAR data
                y_translation_func:
                (int or float)
                                                y coordinate that is used to transform LiDAR data
                center_x_func:
                (float)
                                                Coordinate value of x center.
                                                The x center here was used from the center of
                                                reference DEM without padding
                center_y_func:
                (float)
                                                Coordinate value of x center.
                                                The x center here was used from the center of
                                                reference DEM without padding
                point_classification_func:
                (dictionary)
                                                A dictionary contains point classification of each tile
                point_coordinate_func:
                (dictionary)
                                                A dictionary contains point coordinate (x, y, z) of each tile
                tile_url_list:
                (list)
                                                A list of url
                number_simulation:
                (int)
                                                Ordinal number of simulation
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed_lidar_path = rotated_lidar_path
        transformed = "rotated"
    elif transformation_selection == 't':
        transformed_lidar_path = translated_lidar_path
        transformed = "translated"
    else:
        transformed_lidar_path = combined_lidar_path
        transformed = "combined"

    # Get list of tiles names
    list_tile_paths = glob.glob(fr"{original_lidar_path}\\lidar_{lidar_number}\\Wellington_2013\*.laz")
    length_list = len(list_tile_paths)
    list_tile_names = [pathlib.Path(list_tile_paths[name]).stem for name in range(length_list)]

    # Geometry list of tile index
    tile_index_list = []

    # Transform each tile
    for each_tile in range(length_list):

        # Start timing the whole transformation process -----------------------------------------------------------
        start_transformation_process = time.time()

        # Start timing the transforming process
        start_transforming = time.time()
        # Transform tiles
        rotated_array_func = wrapping_point_rotation(point_coordinate_func[f"{list_tile_names[each_tile]}.laz"],
                                                     angle_func, center_x_func, center_y_func, 0)
        combined_array_func = wrapping_point_translation(rotated_array_func, x_translation_func, y_translation_func)
        # End timing the transforming process
        end_transforming = time.time()

        # Start timing the process of writing LiDAR tiles into las/laz file
        start_writing_lidar = time.time()
        # Write LiDAR tiles into las/laz file
        writing_las_file(transformation_selection, point_classification_func[f"{list_tile_names[each_tile]}.laz"],
                         number_simulation, combined_array_func, list_tile_names[each_tile])
        # End timing the process of writing LiDAR data into las/laz file
        end_writing_lidar = time.time()

        # Start timing the process of writing polygon
        start_writing_polygon = time.time()
        # Write tile index polygon
        tile_polygon = tile_index_polygon(combined_array_func)
        tile_index_list.append(tile_polygon)
        # End timing the process of writing polygon
        end_writing_polygon = time.time()

        # End timing the whole transformation process -------------------------------------------------------------
        end_transformation_process = time.time()

        # PRINT ---------------------------------------------------------------------------------------------------
        # Print separator
        print("*", "Tile", each_tile, "-", list_tile_names[each_tile], "-" * 30)

        # Get the running time of transforming process
        transforming_time = end_transforming - start_transforming

        # Get the running time of writing LiDAR tiles into las/laz file
        writing_lidar_time = end_writing_lidar - start_writing_lidar

        # Get the running time of writing polygon
        writing_polygon = end_writing_polygon - start_writing_polygon

        # Get the running time of whole transformation process
        writing_transformation_process_time = end_transformation_process - start_transformation_process

        # Time list
        time_list = [
            transforming_time / 60,
            writing_lidar_time / 60,
            writing_polygon / 60,
            writing_transformation_process_time / 60
        ]

        # Time title
        time_title = [
            "* Running time of transforming process:",
            "* Running time of writing LiDAR tiles into las/laz file:",
            "* Running time of writing polygon:",
            "* Running time of the whole process:"
        ]

        # Print
        for line in range(len(time_list)):
            print('{0:<60}{1:>10}'.format(time_title[line], time_list[line]))
        print("-" * 30)
        # End printing -------------------------------------------------------------------------------------------

    # Create geopandas dataframe containing tile index information
    tile_info = {
        "Filename": [os.path.basename(list_tile_paths[name]) for name in range(length_list)],
        "URL": tile_url_list
    }
    tile_data_index = gpd.GeoDataFrame(data=tile_info,
                                       geometry=tile_index_list,
                                       crs=2193)

    # Create tiles shapefile
    tile_index_path = pathlib.Path(
        fr"{transformed_lidar_path}\\{transformed}_lidar_{number_simulation}\\Wellington_2013\\Wellington_2013_TileIndex")
    tile_data_index.to_file(tile_index_path)
    shutil.make_archive(base_name=tile_index_path, format='zip', root_dir=tile_index_path)
    shutil.rmtree(tile_index_path)