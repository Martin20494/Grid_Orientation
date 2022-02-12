# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                        # For paths of sub-folders

# Package for paths manipulation
import os                                                   # For manipulating paths and to execute commands
import pathlib                                              # For manipulating the directory path

# Packages for transformation
from transformation import wrapping_point_rotation, \
                           wrapping_point_translation       # For transforming polygon boundaries

# Packages for converting lidar point cloud las/laz file into DEM raster netcdf file
import numpy as np                                          # For all calculation and data array/matrices manipulation
import shapely.geometry                                     # For creating a polygon object
import shutil                                               # For copying file/folder
import xarray                                               # For creating and manipulating an array of raster
import rioxarray                                            # For manipulating spatial data under xarray array format
import json                                                 # For creating json text data format
from geofabrics import processor                            # For executing the conversion of LiDAR into a DEM raster

# Packages for unrotating and untranslating
import rasterio.features                                    # For vectorising features in array
from shapely.geometry import Polygon                        # For creating polygons to change raster values
import geopandas as gpd                                     # For manipulating shape files

# Package for manipulating spatial data
import rasterio                                             # For reading and manipulating spatial data
# ----------------------------------------------------------------------------------------------------------------------


# LIDAR-DERIVED DEMS ###################################################################################################
def dem_raster(transformation_selection,
               resolution_func, chunk_size_func, processor_func,
               number_simulation, padding_func, lidar_dataset_name):
    """ This function is to create a raster file from a las/laz file

    -----------
    References: https://github.com/rosepearson/GeoFabrics
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                resolution_func:
                (int or float)
                                            resolution value in meter
                chunk_size_func:
                (int)
                                            Size value of chunk
                processor_func:
                (int)
                                            Number of processor
                number_simulation:
                (int)
                                            Ordinal number of simulation
                lidar_dataset_name:
                (string)
                                            LiDAR name
                padding_func:
                (list)
                                            A list of x min, x max, y min and y max

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        transformed_lidar_path = rotated_lidar_path
        nc_raster_path_func = rotated_nc_raster_path
        elements_path_func = rotated_elements_path
        transformed = "rotated"
    elif transformation_selection == 't':
        transformed_lidar_path = translated_lidar_path
        nc_raster_path_func = translated_nc_raster_path
        elements_path_func = translated_elements_path
        transformed = "translated"
    else:
        transformed_lidar_path = combined_lidar_path
        nc_raster_path_func = combined_nc_raster_path
        elements_path_func = combined_elements_path
        transformed = "combined"

    # Set crs
    h_crs = 2193
    v_crs = 7839

    # Create folder for storing data materials
    data_dir = pathlib.Path(os.getcwd()) / pathlib.Path(f"{elements_path_func}\\data_{transformed}_{number_simulation}")
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    # Bounding box - Part 1
    x0 = padding_func[0]
    y0 = padding_func[1]
    x1 = padding_func[2]
    y1 = padding_func[3]
    catchment = shapely.geometry.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    catchment = gpd.GeoSeries([catchment])
    catchment = catchment.set_crs(h_crs)

    # Bounding box - Part 2
    catchment_path = pathlib.Path(f"{elements_path_func}\\data_{transformed}_{number_simulation}\\catchment_boundary")
    catchment.to_file(catchment_path)
    shutil.make_archive(base_name=catchment_path, format='zip', root_dir=catchment_path)
    shutil.rmtree(catchment_path)

    # Reference DEM
    grid_dem_x, grid_dem_y = np.meshgrid(np.arange(np.floor(x0), np.ceil(x1), resolution_func),
                                         np.arange(np.floor(y0), np.ceil(y1), resolution_func))
    grid_dem_z = np.full(grid_dem_x.shape, fill_value=0, dtype=np.float64)
    dem = xarray.DataArray(grid_dem_z, coords={'x': grid_dem_x[0], 'y': grid_dem_y[:, 0]}, dims=['y', 'x'],
                           attrs={'scale_factor': 1.0, 'add_offset': 0.0})
    dem.rio.write_crs(h_crs)
    dem.rio.to_raster(f"{elements_path_func}\\data_{transformed}_{number_simulation}\\reference_dem.nc")

    # Design JSON instructions
    instructions = {"instructions":
        {
            "output": {
                "crs": {
                    "horizontal": h_crs,
                    "vertical": v_crs
                },
                "grid_params": {
                    "resolution": resolution_func
                }
            },
            "apis": {
                "open_topography": {
                    f"{lidar_dataset_name}": {
                        "crs": {
                            "horizontal": h_crs,
                            "vertical": v_crs
                        }
                    }
                }
            },
            "processing": {
                "chunk_size": chunk_size_func,
                "number_of_cores": processor_func
            },
            "data_paths": {
                "local_cache": f"{transformed_lidar_path}\\{transformed}_lidar_{number_simulation}",
                "catchment_boundary": f"{elements_path_func}\\data_{transformed}_{number_simulation}\\catchment_boundary.zip",
                "land": f"{elements_path_func}\\data_{transformed}_{number_simulation}\\catchment_boundary.zip",
                "reference_dems": [f"{elements_path_func}\\data_{transformed}_{number_simulation}\\reference_dem.nc"],
                "dense_dem_extents": f"{elements_path_func}\\data_{transformed}_{number_simulation}\\extents.geojson",
                "result_dem": f"{nc_raster_path_func}\\generated_dem_{transformed}_{number_simulation}.nc"
            },
            "general": {
                "set_dem_shoreline": True,
                "drop_offshore_lidar": False,
                "lidar_classifications_to_keep": [2],
                "interpolation_missing_values": True
            }
        }
    }

    # Create json instructions
    with open(f"{elements_path_func}\\data_{transformed}_{number_simulation}\\instructions.json", "w") as instruction:
        json.dump(instructions, instruction)

    # Create DEM raster
    runner = processor.DemGenerator(instructions)
    runner.run()


# END LIDAR-DERIVED DEMS ###############################################################################################


# POLYGON BOUNDARIES ###################################################################################################
def polygon_generation(polygon_boundary, transformation_selection, number_simulation,
                       angle_func, x_translation_func, y_translation_func,
                       center_x_func, center_y_func,
                       filename_output_polygon):
    """This function is to create polygon layer to change the pixel values

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                polygon_boundary:
                                                A list of coordinates of
                transformation_selection:
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                number_simulation:
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
                filename_output_polygon:
                (string)
                                                Output file name
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        tiff_raster_path_func = rotated_tiff_raster_path
        transformed = "rotated"
    elif transformation_selection == 't':
        tiff_raster_path_func = translated_tiff_raster_path
        transformed = "translated"
    else:
        tiff_raster_path_func = combined_tiff_raster_path
        transformed = "combined"

    # Create storing folder to store the file
    polygon_dir = pathlib.Path(os.getcwd()) / pathlib.Path(
        fr"{tiff_raster_path_func}\\{transformed}_{number_simulation}")
    if not os.path.exists(polygon_dir):
        os.mkdir(polygon_dir)

    # Convert the list into array to do transformation
    boundary_coordinates = np.array(polygon_boundary).astype('float64')

    # Transform
    rotated_boundary = wrapping_point_rotation(boundary_coordinates, angle_func, center_x_func, center_y_func, 0)
    translated_boundary = wrapping_point_translation(rotated_boundary, x_translation_func, y_translation_func)
    transformed_boundary = translated_boundary

    # Convert into shapely geometry Polygon
    shapely_polygon = Polygon(transformed_boundary)

    # Create geopandas dataframe
    polygon_dataframe = gpd.GeoDataFrame(geometry=[shapely_polygon],
                                         crs=2193)

    # Write into shape file
    polygon_dataframe.to_file(fr"{filename_output_polygon}.shp")


def polygon_boundaries(transformation_selection, number_simulation,
                       resolution_func,
                       angle_func, x_translation_func, y_translation_func,
                       center_x_func, center_y_func):
    """This function is to create polygon boundaries by using polygon_generation() function.
       These polygons will be used to change the pixel values inside or outside.

    -----------
    References: https://gis.stackexchange.com/questions/414194/changing-raster-pixel-values-outside-of-polygon-box-using-rioxarray
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
                resolution_func:
                (int or float)
                                                Resolution value in meter
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
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up the path for transformation_selection
    if transformation_selection == 'r':
        nc_raster_path_func = rotated_nc_raster_path
        tiff_raster_path_func = rotated_tiff_raster_path
        transformed = "rotated"
    elif transformation_selection == 't':
        nc_raster_path_func = translated_nc_raster_path
        tiff_raster_path_func = translated_tiff_raster_path
        transformed = "translated"
    else:
        nc_raster_path_func = combined_nc_raster_path
        tiff_raster_path_func = combined_tiff_raster_path
        transformed = "combined"

    # Define the path for output
    output_path = fr"{tiff_raster_path_func}\\{transformed}_{number_simulation}"

    # PADDING #####################################################
    # Get extent information of original no-padding raster
    raster_no_padding = rasterio.open(fr"{nc_raster_path_func}\\generated_dem_no_padding.nc")
    no_padding_xmin, no_padding_ymin, no_padding_xmax, no_padding_ymax = raster_no_padding.bounds

    # Create polygon to remove PADDING values
    # Set up polygon coordinates into a list
    polygon_padding = [(no_padding_xmax, no_padding_ymax),
                       (no_padding_xmax, no_padding_ymin),
                       (no_padding_xmin, no_padding_ymin),
                       (no_padding_xmin, no_padding_ymax),
                       (no_padding_xmax, no_padding_ymax)]

    # Create polygon padding layer under shape file
    polygon_generation(polygon_padding,
                       transformation_selection,
                       number_simulation,
                       angle_func, x_translation_func, y_translation_func,
                       center_x_func, center_y_func,
                       fr"{output_path}\\polygon_padding_{number_simulation}")

    # BOUNDARY 1 ---------------------------------------------------------------------------------------------
    # Get extent information for first boundary

    # # Boundary for test case 1
    # boundary_1_ymin = no_padding_ymin + resolution_func*2
    # boundary_1_xmax = no_padding_xmax
    # boundary_1_ymax = boundary_1_ymin + 10*7
    # boundary_1_xmin = no_padding_xmax - resolution_func

    # # Boundary for test case 2
    # boundary_1_ymin = no_padding_ymin + resolution_func*1
    # boundary_1_xmax = no_padding_xmax
    # boundary_1_ymax = boundary_1_ymin + 10*3
    # boundary_1_xmin = no_padding_xmax - resolution_func

    # # Boundary for large area
    # boundary_1_ymin = no_padding_ymin + 10*59
    # boundary_1_xmax = no_padding_xmax
    # boundary_1_ymax = boundary_1_ymin + 10*10
    # boundary_1_xmin = no_padding_xmax - 10

    # Boundary for small area
    boundary_1_ymin = no_padding_ymin + resolution_func * (10*12/resolution_func)
    boundary_1_xmax = no_padding_xmax
    boundary_1_ymax = boundary_1_ymin + resolution_func * (10*22/resolution_func)
    boundary_1_xmin = no_padding_xmax - resolution_func

    # # Boundary for Buller
    # boundary_1_ymin = no_padding_ymin
    # boundary_1_ymax = no_padding_ymin + resolution_func
    # boundary_1_xmin = no_padding_xmin + resolution_func * (750/resolution_func)
    # boundary_1_xmax = boundary_1_xmin + resolution_func * (560/resolution_func)

    # Write into a list
    polygon_boundary_1 = [(boundary_1_xmax, boundary_1_ymax),
                          (boundary_1_xmax, boundary_1_ymin),
                          (boundary_1_xmin, boundary_1_ymin),
                          (boundary_1_xmin, boundary_1_ymax),
                          (boundary_1_xmax, boundary_1_ymax)]

    # Create polygon boundary 1 under shape file
    polygon_generation(polygon_boundary_1,
                       transformation_selection,
                       number_simulation,
                       angle_func, x_translation_func, y_translation_func,
                       center_x_func, center_y_func,
                       fr"{output_path}\\polygon_boundary_1_{number_simulation}")
    # ---------------------------------------------------------------------------------------------------------

    # BOUNDARY 2 ----------------------------------------------------------------------------------------------
    # Get extent information from raster including padding
    raster_padding = rasterio.open(fr"{nc_raster_path_func}\\generated_dem_reference.nc")
    padding_xmin, padding_ymin, padding_xmax, padding_ymax = raster_padding.bounds

    # Get coordinates of boundary 2
    boundary_2_ymax = no_padding_ymin
    boundary_2_xmax = padding_xmax
    boundary_2_ymin = padding_ymin
    boundary_2_xmin = padding_xmin

    # Write into a list
    polygon_boundary_2 = [(boundary_2_xmax, boundary_2_ymax),
                          (boundary_2_xmax, boundary_2_ymin),
                          (boundary_2_xmin, boundary_2_ymin),
                          (boundary_2_xmin, boundary_2_ymax),
                          (boundary_2_xmax, boundary_2_ymax)]

    # Create polygon boundary 2 under shape file
    polygon_generation(polygon_boundary_2,
                       transformation_selection,
                       number_simulation,
                       angle_func, x_translation_func, y_translation_func,
                       center_x_func, center_y_func,
                       fr"{output_path}\\polygon_boundary_2_{number_simulation}")
    # ---------------------------------------------------------------------------------------------------------


# END POLYGON BOUNDARIES ###############################################################################################


# VALUE CHANGES ########################################################################################################
def value_change(shapefile_func, file_need_changing_func, value_func, inside=True):
    """This function is to change pixel values inside or outside polygons

    -----------
    References: https://corteva.github.io/rioxarray/html/rioxarray.html
                https://corteva.github.io/rioxarray/stable/examples/convert_to_raster.html
                https://gis.stackexchange.com/questions/414194/changing-raster-pixel-values-outside-of-polygon-box-using-rioxarray
                https://automating-gis-processes.github.io/CSC/notebooks/L2/geopandas-basics.html
                https://corteva.github.io/rioxarray/stable/getting_started/nodata_management.html
                https://gdal.org/programs/gdal_translate.html
                https://gis.stackexchange.com/questions/390438/replacing-nodata-values-by-a-constant-using-gdal-cli-tools
    -----------

    -----------
    Arguments:
                shapefile_func:
                (polygon)
                                            Polygon boundaries
                file_need_changing_func:
                (string)
                                            File contains values that need changing
                value_func:
                (int or float)
                                            Values used to replace
                inside:
                (boolean)
                                            If True, change values inside, else, change values outside
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up value changing command
    if inside:
        # Change values inside polygons
        inside_command = fr"gdal_rasterize -burn {value_func} {shapefile_func} {file_need_changing_func}"
        os.system(inside_command)
    else:
        # Change values outside polygons
        outside_command = fr"gdal_rasterize -i -burn {value_func} {shapefile_func} {file_need_changing_func}"
        os.system(outside_command)


def value_change_execution(transformation_selection, number_simulation):
    """ This function is to apply value_change() function on a certain raster

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
        nc_raster_path_func = rotated_nc_raster_path
        tiff_raster_path_func = rotated_tiff_raster_path
        transformed = "rotated"
    elif transformation_selection == 't':
        nc_raster_path_func = translated_nc_raster_path
        tiff_raster_path_func = translated_tiff_raster_path
        transformed = "translated"
    else:
        nc_raster_path_func = combined_nc_raster_path
        tiff_raster_path_func = combined_tiff_raster_path
        transformed = "combined"

    # Convert NetCDF file into GeoTiff file
    raster_nc = rioxarray.open_rasterio(fr"{nc_raster_path_func}\\generated_dem_{transformed}_{number_simulation}.nc")
    raster_nc.rio.to_raster(
        fr"{tiff_raster_path_func}\\{transformed}_{number_simulation}\\generated_dem_{transformed}_{number_simulation}.tif")

    # Get polygon shape_file and file that has values need changing
    file_path = fr"{tiff_raster_path_func}\\{transformed}_{number_simulation}"
    shapefile_padding_dir = fr"{file_path}\\polygon_padding_{number_simulation}.shp"
    shapefile_boundary_1_dir = fr"{file_path}\\polygon_boundary_1_{number_simulation}.shp"
    shapefile_boundary_2_dir = fr"{file_path}\\polygon_boundary_2_{number_simulation}.shp"
    file_need_changing_dir = fr"{file_path}\\generated_dem_{transformed}_{number_simulation}.tif"

    # Changing the value
    # Changing padding values into -999
    value_change(shapefile_padding_dir, file_need_changing_dir, -999, False)

    # Changing boundary 1 values into 100
    value_change(shapefile_boundary_1_dir, file_need_changing_dir, 100, True)

    # Changing boundary 2 values into 100
    value_change(shapefile_boundary_2_dir, file_need_changing_dir, 100, True)


# END VALUES CHANGE ####################################################################################################


# ASCII FILE CONVERSION ################################################################################################
def convert_to_asc(transformation_selection, number_simulation):
    """This function is to convert GeoTiff file into ASCII file

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
        tiff_raster_path_func = rotated_tiff_raster_path
        asc_raster_path_func = rotated_asc_raster_path
        transformed = "rotated"
    elif transformation_selection == 't':
        tiff_raster_path_func = translated_tiff_raster_path
        asc_raster_path_func = translated_asc_raster_path
        transformed = "translated"
    else:
        tiff_raster_path_func = combined_tiff_raster_path
        asc_raster_path_func = combined_asc_raster_path
        transformed = "combined"

    # Convert GeoTiff file into ASCII file
    tiff_raster_func = rioxarray.open_rasterio(
        fr"{tiff_raster_path_func}\\{transformed}_{number_simulation}\\generated_dem_{transformed}_{number_simulation}.tif")
    nodata_tiff_raster_func = tiff_raster_func.rio.write_nodata(-999, inplace=True)
    crs_tiff_raster_func = nodata_tiff_raster_func.rio.write_crs(2193)
    crs_tiff_raster_func.rio.to_raster(fr"{asc_raster_path_func}\\generated_dem_{transformed}_{number_simulation}.asc")

# END ASCII FILE CONVERSION ############################################################################################
