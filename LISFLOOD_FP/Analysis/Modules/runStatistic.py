# Prepare packages -----------------------------------------------------------------------------------------------------

from statisticalCalculation import statistic_calculation, \
                                   area_calculation, \
                                   building_calculation                 # For statistical calculation

from fileWriting import raster_generation                               # For raster generation


# ----------------------------------------------------------------------------------------------------------------------

def calculation_dict(transformation_selection, dataset_func,
                     resolution_func, filter_rate, building_path, dif_path=None, raster=True,
                     rectangle=True, switch=False):
    """

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
                dataset_func:
                (pandas dataframe)
                                                A full dataset including necessary information
                resolution_func:
                (int)
                                                Resolution in meter
                filter_rate:
                (float or int)
                                                The rate at which the depth values will be ignored
                building_path:
                (string)
                                                Path of file containing building polygons
                dif_path:
                (string)
                                                None if dif_path is not specified
                                                Different path for different file version
                raster:
                (boolean)
                                                To decide if raster should be generated
                rectangle:
                (boolean)
                                                To specify that if the data shape is rectangle or not.
                                                True is rectangle (default)
                                                False is not rectangle
                switch:
                (boolean)
                                                To switch the row and column of array
                                                True is switching
                                                False is not switching (default)
    -----------

    -----------
    Returns:
                calculation_dict_set:
                (dictionary)
                                                A dictionary contains all statistical results
    -----------

    """

    statistic_method = ["mean", "sd", "cv", "cell"]

    calculation_dict_set = {}

    for each_method in statistic_method:
        # Statistical calculation
        statistic_result = statistic_calculation(dataset_func, each_method, filter_rate)

        # Store in dictionary
        calculation_dict_set[each_method] = statistic_result

        # Write into raster
        if raster:
            raster_generation(
                transformation_selection,
                statistic_result['x'],
                statistic_result['y'],
                statistic_result[each_method],
                each_method,
                None,
                rectangle,
                switch
            )
        else:
            pass

    # Add area result
    calculation_dict_set['area'] = area_calculation(dataset_func, resolution_func, filter_rate)

    # Add building result
    calculation_dict_set['building'] = building_calculation(transformation_selection,
                                                            dataset_func, building_path, dif_path)

    return calculation_dict_set
