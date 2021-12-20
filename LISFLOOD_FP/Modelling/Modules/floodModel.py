# Prepare packages -----------------------------------------------------------------------------------------------------
# Packages for path manipulation
from folder import *                                       # For paths of sub-folders
import os                                                  # For manipulating the directory path and to execute commands
import pathlib                                             # For manipulating the directory path
import shutil                                              # For copying file/folder

# Packages for transformation
import numpy as np                                          # For all calculation and data array/matrices manipulation
from transformation import wrapping_point_rotation, \
                           wrapping_point_translation       # For transforming polygon boundaries

# ----------------------------------------------------------------------------------------------------------------------


def parameter_files(transformation_selection, number_simulation,
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
    #     point_source = np.array([(1772787.2, 5472314.7)]).astype('float64')

    point_source = np.array([(1773016.362, 5472061.027)]).astype('float64')

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
    discharge = [(25, 0),
                 (40, 300),
                 (40, 600),
                 (45, 1200),
                 (50, 1800),
                 (55, 2400),
                 (60, 3000),
                 (65, 3600)]
    discharge_array = np.array(discharge)
    row_number = discharge_array.shape[0]

    # Write into text file bci format
    with open(fr"{param_dir}\\{transformed}_{number_simulation}.bdy", "w") as discharge:
        note = 'Waikanae - Predict in 1 hour\n'
        reference = note + 'Waikanae\n'
        column_name = reference + '{0:<8}seconds\n'.format(row_number)
        discharge.write(column_name)
        for line in range(row_number):
            data_discharge = discharge_array[line]
            text_river_discharge = '{0[0]:<8}{0[1]}\n'.format(data_discharge)
            discharge.write(text_river_discharge)

    # Create output path
    # Create folder to store the parameter files
    output_dir = pathlib.Path(os.getcwd()) / pathlib.Path(fr"{output_path}\\{transformed}_{number_simulation}")
    if not os.path.exists(param_dir):
        os.mkdir(param_dir)

    # PAR FILE ---------------------------------------------------------------
    # Construct parameter files
    parameters_list = [('resroot', 'out'),
                       ('saveint', 300),
                       ('massint', 100),
                       ('sim_time', 3600),
                       ('initial_tstep', 2),
                       ('bcifile', fr"{param_dir}\\{transformed}_{number_simulation}.bci"),
                       ('bdyfile', fr"{param_dir}\\{transformed}_{number_simulation}.bdy"),
                       ('DEMFile', fr"{asc_raster_path_func}\\generated_dem_{transformed}_{number_simulation}.asc"),
                       ('fpfric', 0.06),
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


def run_LISFLOOD(transformation_selection, number_simulation,
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
        parameter_files('r', number_simulation,
                        angle_func, x_translation_func, y_translation_func,
                        center_x_func, center_y_func)
        os.chdir(r'P:\\Courses\\PhD\\LISFLOOD_FP')
        os.system("lisflood_v8.exe -v Waikanae_LISFLOOD_acceleration.par")
    elif transformation_selection == "t":
        parameter_files('t', number_simulation,
                        angle_func, x_translation_func, y_translation_func,
                        center_x_func, center_y_func)
        os.chdir(r'P:\\Courses\\PhD\\LISFLOOD_FP')
        os.system("lisflood_v8.exe -v Waikanae_LISFLOOD_acceleration.par")
    else:
        parameter_files('c', number_simulation,
                        angle_func, x_translation_func, y_translation_func,
                        center_x_func, center_y_func)
        os.chdir(r'P:\\Courses\\PhD\\LISFLOOD_FP')
        os.system("lisflood_v8.exe -v Waikanae_LISFLOOD_acceleration.par")

    # Change the path back to origin
    os.chdir(fr"{header}\\")