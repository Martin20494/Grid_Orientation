# Prepare packages -----------------------------------------------------------------------------------------------------
# Packages for path manipulation
from folder import *                                       # For paths of sub-folders
import os                                                  # For manipulating the directory path and to execute commands
import pathlib                                             # For manipulating the directory path
import shutil                                              # For copying file/folder

# Packages for transformation
import numpy as np                                         # For all calculation and data array/matrices manipulation
from transformation import wrapping_point_rotation, \
                           wrapping_point_translation      # For transforming polygon boundaries

# Packages for model
from hydrograph import discharge_data                      # For generating discharge following to resolution
# ----------------------------------------------------------------------------------------------------------------------


def parameter_files(transformation_selection,
                    discharge_file, time_file, resolution_func,
                    number_simulation,
                    angle_func, x_translation_func, y_translation_func,
                    center_x_func, center_y_func):
    """This function is to create parameters files for LISFLOOD-FP model

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
                discharge_file:
                (string)
                                                Path to the discharge file

                time_file:
                (string)
                                                Path to the time file
                resolution_func:
                (int or float)
                                                resolution value in meter
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
        asc_raster_path_func = rotated_asc_raster_path
        output_path = rotated_FPoutput_path
        transformed = "rotated"
    elif transformation_selection == 't':
        asc_raster_path_func = translated_asc_raster_path
        output_path = translated_FPoutput_path
        transformed = "translated"
    else:
        asc_raster_path_func = combined_asc_raster_path
        output_path = combined_FPoutput_path
        transformed = "combined"

    # BCI FILE --------------------------------------------------------------
    # Get the point source

    # point_source = np.array([(1483758.308, 5373080.117)]).astype('float64') # Buller River

    # point_source = np.array([(1773704.959, 5472333.515)]).astype('float64') # larger area

    point_source = np.array([(1773016.362, 5472061.027)]).astype('float64') # smaller area

    # point_source = np.array([(1769951.132, 5472908.887)]).astype('float64')  # test case 1

    # point_source = np.array([(1769795.781, 5472903.947)]).astype('float64') # test case 2

    # Transform the point source
    rotated_point_source = wrapping_point_rotation(point_source, angle_func, center_x_func, center_y_func, 0)
    translated_point_source = wrapping_point_translation(rotated_point_source, x_translation_func, y_translation_func)
    transformed_point_source = translated_point_source

    # Construct floodplain boundary information
    floodplain_boundary = ['P', str(transformed_point_source[0, 0]), str(transformed_point_source[0, 1]), 'QVAR',
                           'Waikanae']

    # Create folder to store the parameter files
    param_dir = pathlib.Path(os.getcwd()) / pathlib.Path(
        fr"P:\\Courses\\PhD\\LISFLOOD_FP\\{transformed}_param\\{transformed}_{number_simulation}")
    if not os.path.exists(param_dir):
        os.mkdir(param_dir)

    # Write into text file bci format
    with open(fr"{param_dir}\\{transformed}_{number_simulation}.bci", "w") as fp_boundary:
        text_fp_boundary = '{0[0]:<5}{0[1]:<20}{0[2]:<20}{0[3]:<7}{0[4]:<5}'.format(floodplain_boundary)
        fp_boundary.write(text_fp_boundary)

    # BDY FILE --------------------------------------------------------------
    # Construct river discharge
    discharge_array = discharge_data(discharge_file, time_file, resolution_func)
    row_number = discharge_array.shape[0]

    # Write into text file bci format
    with open(fr"{param_dir}\\{transformed}_{number_simulation}.bdy", "w") as discharge:
        note = 'Waikanae - Predict in 1 hour\n'
        reference = note + 'Waikanae\n'
        column_name = reference + '{0:<20}seconds\n'.format(row_number)
        discharge.write(column_name)
        for line in range(row_number):
            data_discharge = discharge_array[line]
            text_river_discharge = '{0[0]:<20}{0[1]}\n'.format(data_discharge)
            discharge.write(text_river_discharge)

    # Create output path
    # Create folder to store the parameter files
    output_dir = pathlib.Path(os.getcwd()) / pathlib.Path(fr"{output_path}\\{transformed}_{number_simulation}")
    if not os.path.exists(param_dir):
        os.mkdir(param_dir)

    # PAR FILE ---------------------------------------------------------------
    # Construct parameter files
    parameters_list = [('resroot', 'out'),
                       ('saveint', 100),       #for normal saveint = 200 if sim_time = 7200. Real event, saveint = 100
                       ('massint', 100),       #for normal massint = 100 if sim_time = 7200. Real event, massint = 100
                       ('sim_time', 4800),     #for normal sim_time = 7200. Real event, sim_time = 4800
                       ('initial_tstep', 2),
                       ('bcifile', fr"{param_dir}\\{transformed}_{number_simulation}.bci"),
                       ('bdyfile', fr"{param_dir}\\{transformed}_{number_simulation}.bdy"),
                       ('DEMFile', fr"{asc_raster_path_func}\\generated_dem_{transformed}_{number_simulation}.asc"),
                       ('fpfric', 0.03),
                       ('dirroot', fr"{output_dir}")]
    parameters_array = np.array(parameters_list)

    # Write into text file par format
    with open(fr"P:\\Courses\\PhD\\LISFLOOD_FP\\Waikanae_LISFLOOD_acceleration.par", 'w') as parameters:
        for each_parameter in range(parameters_array.shape[0]):
            data_parameter = parameters_array[each_parameter]
            text_parameter = '{0[0]:<20}{0[1]}\n'.format(data_parameter)
            parameters.write(text_parameter)
        parameters.write('acceleration\ndrain_data\n\n')

    # Copy the text files for later checking
    shutil.copy2(fr"P:\\Courses\\PhD\\LISFLOOD_FP\\Waikanae_LISFLOOD_acceleration.par",
                 fr"{param_dir}\\{transformed}_{number_simulation}.par")


def run_LISFLOOD(transformation_selection,
                 discharge_file, time_file, resolution_func,
                 number_simulation,
                 angle_func, x_translation_func, y_translation_func,
                 center_x_func, center_y_func):
    """This function is to run LISFLOOD-FP model

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
                discharge_file:
                (string)
                                            Path to the discharge file

                time_file:
                (string)
                                            Path to the time file
                resolution_func:
                (int or float)
                                            resolution value in meter
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
    if transformation_selection == "r":
        parameter_files('r',
                        discharge_file, time_file, resolution_func,
                        number_simulation,
                        angle_func, x_translation_func, y_translation_func,
                        center_x_func, center_y_func)
        os.chdir(r'P:\\Courses\\PhD\\LISFLOOD_FP')
        os.system("lisflood_v8.exe -v Waikanae_LISFLOOD_acceleration.par")
    elif transformation_selection == "t":
        parameter_files('t',
                        discharge_file, time_file, resolution_func,
                        number_simulation,
                        angle_func, x_translation_func, y_translation_func,
                        center_x_func, center_y_func)
        os.chdir(r'P:\\Courses\\PhD\\LISFLOOD_FP')
        os.system("lisflood_v8.exe -v Waikanae_LISFLOOD_acceleration.par")
    else:
        parameter_files('c',
                        discharge_file, time_file, resolution_func,
                        number_simulation,
                        angle_func, x_translation_func, y_translation_func,
                        center_x_func, center_y_func)
        os.chdir(r'P:\\Courses\\PhD\\LISFLOOD_FP')
        os.system("lisflood_v8.exe -v Waikanae_LISFLOOD_acceleration.par")

    # Change the path back to origin
    os.chdir(fr"{header}\\")