# Prepare packages -----------------------------------------------------------------------------------------------------

from statisticalCalculation import statistic_calculation, \
                                   area_calculation                     # For statistical calculation

from fileWriting import raster_generation                               # For raster generation


# ----------------------------------------------------------------------------------------------------------------------

def calculation_dict(transformation_selection, dataset_func, resolution_func, filter_rate):
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
        raster_generation(
            transformation_selection,
            statistic_result['x'],
            statistic_result['y'],
            statistic_result[each_method],
            each_method
        )

    # Add area result
    calculation_dict_set['area'] = area_calculation(dataset_func, resolution_func, filter_rate)

    return calculation_dict_set
