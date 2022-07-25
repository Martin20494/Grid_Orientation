# Prepare packages -----------------------------------------------------------------------------------------------------
from transformation import get_dictionary,\
                           transformed_tiles_generation            # For generating LiDAR simulations of transformation
from demRaster import dem_raster, \
                      polygon_boundaries, \
                      value_change_execution, \
                      convert_to_asc                               # For creating DEMs
from floodModel import run_LISFLOOD                                # For flood model
from flowdepthExtraction import flowdepth_extraction               # For extract flowdepth
from untransformation import polygons_untransformation             # For untransformation
import time                                                        # For timing steps
# ----------------------------------------------------------------------------------------------------------------------


def run_simulation(set_of_simulation,
                   transformation_selection, lidar_number,
                   center_x_func, center_y_func,
                   tile_url_list, lidar_dataset_name):
    """This function is to run transformation process with given information simulations

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                set_of_simulation:
                (array)
                                                An 3D array of simulations (angle, x, y)
                transformation_selection:
                (string)
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                lidar_number:
                (int)
                                                Ordinal number of lidar folder
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
                tile_url_list:
                (list)
                                                A list of url
                lidar_dataset_name:
                (string)
                                                LiDAR name
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Get dictionaries of classification and coordinates
    point_classification = get_dictionary(lidar_number, "classification", lidar_dataset_name)
    point_coordinate = get_dictionary(lidar_number, "coordinate", lidar_dataset_name)

    # Iterate each simulation
    for order_number in set_of_simulation:
        # Get values
        angle_val = order_number[0]
        x_val = order_number[1]
        y_val = order_number[2]

        # Set the name of each transformation
        combination_number = f"angle_{angle_val}_x_{x_val}_y_{y_val}"

        # Print the header of the first running
        print("{0} angle: {1} and x: {2} and y: {3} {0}".format("-" * 30, angle_val, x_val, y_val))

        # Start timing transforming process -------------------------------------------------------------
        start_transform_tiles = time.time()
        # Transform tiles
        transformed_tiles_generation(
            transformation_selection, lidar_number,
            angle_val, x_val, y_val,
            center_x_func, center_y_func,
            tile_url_list,
            point_classification, point_coordinate,
            combination_number, lidar_dataset_name
        )
        # End timing generating DEM process ------------------------------------------------------------
        end_transform_tiles = time.time()

        # PRINT -----------------------------------------------------------------------------------------
        # Get the running time of transforming tiles process
        transform_tile_time = end_transform_tiles - start_transform_tiles

        # Convert second to minute
        transform_tile_minute = transform_tile_time / 60

        # Get the name
        transform_tile_name = "* Running time of transforming tiles:"

        # Print
        print('{0:<60}{1:>10}'.format(transform_tile_name, transform_tile_minute))

        # Print separator
        print("-" * 90)
        print("\n")


def run_dem(set_of_simulation,
            transformation_selection,
            resolution_func, chunk_size_func, processor_func,
            dem_boundary, lidar_dataset_name):
    """This function is to run dem generating process

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                set_of_simulation:
                (array)
                                                An 3D array of simulations (angle, x, y)
                transformation_selection:
                (string)
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                resolution_func:
                (int)
                                                resolution value in meter
                chunk_size_func:
                (int)
                                                Size value of chunk
                processor_func:
                (int)
                                                Number of processor
                dem_boundary:
                (list)
                                                A list of boundary (xmin, ymin, xmax, ymax)
                lidar_dataset_name:
                (string)
                                                LiDAR dataset name
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Iterate each simulation
    for order_number in set_of_simulation:
        # Get values
        angle_val = order_number[0]
        x_val = order_number[1]
        y_val = order_number[2]

        # Set the name of each transformation
        combination_number = f"angle_{angle_val}_x_{x_val}_y_{y_val}"

        # Print the header of the first running
        print("{0} angle: {1} and x: {2} and y: {3} {0}".format("-" * 30, angle_val, x_val, y_val))

        # Start timing generating DEM process -----------------------------------------------------------
        start_dem = time.time()
        # Generate DEM
        dem_raster(transformation_selection,
                   resolution_func, chunk_size_func, processor_func,
                   combination_number, dem_boundary, lidar_dataset_name)
        # Start timing generating DEM process -----------------------------------------------------------
        end_dem = time.time()

        # PRINT -----------------------------------------------------------------------------------------
        # Get the running time
        dem_time = end_dem - start_dem

        # Convert second to minute
        dem_minute = dem_time / 60

        # Get the name
        dem_name = "* Running time of dem generation:"

        # Print
        print('{0:<60}{1:>10}'.format(dem_name, dem_minute))

        # Print separator
        print("-" * 90)
        print("\n")


def run_changing_value(set_of_simulation,
                       transformation_selection,
                       resolution_func,
                       center_x_func, center_y_func):
    """This function is to run changing values process

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                set_of_simulation:
                (array)
                                                An 3D array of simulations (angle, x, y)
                transformation_selection:
                (string)
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                resolution_func:
                (int or float)
                                                Resolution value in meter
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
    # Iterate each simulation
    for order_number in set_of_simulation:
        # Get values
        angle_val = order_number[0]
        x_val = order_number[1]
        y_val = order_number[2]

        # Set the name of each transformation
        combination_number = f"angle_{angle_val}_x_{x_val}_y_{y_val}"

        # Print the header of the first running
        print("{0} angle: {1} and x: {2} and y: {3} {0}".format("-" * 30, angle_val, x_val, y_val))

        # Start timing the whole process ---------------------------------------------------------------
        start_changing_values = time.time()

        # Create polygon boundaries
        polygon_boundaries(transformation_selection, combination_number,
                           resolution_func,
                           angle_val, x_val, y_val,
                           center_x_func, center_y_func)

        # Changing pixels values
        value_change_execution(transformation_selection, combination_number)

        # Convert to ASCII
        convert_to_asc(transformation_selection, combination_number)

        # End timing the whole process ------------------------------------------------------------------
        end_changing_values = time.time()

        # PRINT -----------------------------------------------------------------------------------------
        # Get the running time
        changing_values_time = end_changing_values - start_changing_values

        # Convert second to minute
        changing_values_minute = changing_values_time / 60

        # Get the name
        changing_values_name = "* Running time of changing values process:"

        # Print
        print('{0:<60}{1:>10}'.format(changing_values_name, changing_values_minute))

        # Print separator
        print("-" * 90)
        print("\n")


def run_flood_model(set_of_simulation,
                    transformation_selection,
                    discharge_file, time_file, resolution_func,
                    center_x_func, center_y_func):
    """This function is to run flood model

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                set_of_simulation:
                (array)
                                                An 3D array of simulations (angle, x, y)
                transformation_selection:
                (string)
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                discharge_file:
                (string)
                                                Path to the discharge file

                time_file:
                (string)
                                                Path to the time file
                resolution_func:
                (int or float)
                                                resolution value in meter
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
    # Iterate each simulation
    for order_number in set_of_simulation:
        # Get values
        angle_val = order_number[0]
        x_val = order_number[1]
        y_val = order_number[2]

        # Set the name of each transformation
        combination_number = f"angle_{angle_val}_x_{x_val}_y_{y_val}"

        # Print the header of the first running
        print("{0} angle: {1} and x: {2} and y: {3} {0}".format("-" * 30, angle_val, x_val, y_val))

        # Start timing the whole model running process --------------------------------------------
        start_model = time.time()

        # Run LISFLOOD-FP
        run_LISFLOOD(transformation_selection,
                     discharge_file, time_file, resolution_func,
                     combination_number,
                     angle_val, x_val, y_val,
                     center_x_func, center_y_func)

        # End timing the whole model running process ----------------------------------------------
        end_model = time.time()

        # PRINT -----------------------------------------------------------------------------------------
        # Get the running time
        model_time = end_model - start_model

        # Convert second to minute
        model_minute = model_time / 60

        # Get the name
        model_name = "* Running time of running flood model:"

        # Print
        print('{0:<60}{1:>10}'.format(model_name, model_minute))

        # Print separator
        print("-" * 90)
        print("\n")


def run_flowdepth_extraction(set_of_simulation,
                             transformation_selection,
                             time_extract_func):
    """This function is to run flowdepth extraction process

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                set_of_simulation:
                (array)
                                                An 3D array of simulations (angle, x, y)
                transformation_selection:
                (string)
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                time_extract_func:
                (int)
                                                Predicted time for running flood model
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Iterate each simulation
    for order_number in set_of_simulation:
        # Get values
        angle_val = order_number[0]
        x_val = order_number[1]
        y_val = order_number[2]

        # Set the name of each transformation
        combination_number = f"angle_{angle_val}_x_{x_val}_y_{y_val}"

        # Print the header of the first running
        print("{0} angle: {1} and x: {2} and y: {3} {0}".format("-" * 30, angle_val, x_val, y_val))

        # Start timing the whole extracting process -----------------------------------------------------
        start_extracting = time.time()

        # Extract flowdepth
        flowdepth_extraction(transformation_selection, combination_number, time_extract_func)

        # End timing the whole extracting process -------------------------------------------------------
        end_extracting = time.time()

        # PRINT -----------------------------------------------------------------------------------------
        # Get the running time
        extracting_time = end_extracting - start_extracting

        # Convert second to minute
        extracting_minute = extracting_time / 60

        # Get the name
        extracting_name = "* Running time of extracting flowdepth:"

        # Print
        print('{0:<60}{1:>10}'.format(extracting_name, extracting_minute))

        # Print separator
        print("-" * 90)
        print("\n")


def run_untransformation(set_of_simulation,
                         transformation_selection,
                         time_extract_func):
    """This function is to run un-transformation process

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                set_of_simulation:
                (array)
                                                An 3D array of simulations (angle, x, y)
                transformation_selection:
                (string)
                                                "r" means rotation
                                                "t" means translation
                                                "c" means combination
                time_extract_func:
                (int)
                                                Predicted time for running flood model
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Iterate each simulation
    for order_number in set_of_simulation:
        # Get values
        angle_val = order_number[0]
        x_val = order_number[1]
        y_val = order_number[2]

        # Set the name of each transformation
        combination_number = f"angle_{angle_val}_x_{x_val}_y_{y_val}"

        # Print the header of the first running
        print("{0} angle: {1} and x: {2} and y: {3} {0}".format("-" * 30, angle_val, x_val, y_val))

        # Translated values for raster
        x_raster = x_val / 10 * (-1)
        y_raster = y_val / 10

        # Start timing the whole un-combination process --------------------------------------------------
        start_untransform = time.time()

        # Un-combination
        polygons_untransformation(transformation_selection, combination_number,
                                  angle_val, x_raster, y_raster,
                                  time_extract_func)

        # End timing the whole un-rotation process -------------------------------------------------------
        end_untransform = time.time()

        # PRINT -----------------------------------------------------------------------------------------
        # Get the running time
        untransform_time = end_untransform - start_untransform

        # Convert second to minute
        untransform_minute = untransform_time / 60

        # Get the name
        untransform_name = "* Running time of untransformation:"

        # Print
        print('{0:<60}{1:>10}'.format(untransform_name, untransform_minute))

        # Print separator
        print("-" * 90)
        print("\n")