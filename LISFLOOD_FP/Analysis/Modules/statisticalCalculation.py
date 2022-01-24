# Prepare packages -----------------------------------------------------------------------------------------------------
import pandas as pd                                                  # For manipulating dataframe
# ----------------------------------------------------------------------------------------------------------------------

def statistic_calculation(dataset_func, calculation_option, filter_rate_func):
    """This function is to do some statistical calculation

    -----------
    References: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.std.html
                https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.mean.html
    -----------

    -----------
    Arguments:
                dataset_func:
                (pandas dataframe)
                                        A full dataset including unnecessary information
                calculation_option:
                (string)
                                        "mean" means to calculate mean
                                        "sd" means to calculate standard deviation
                                        "cv" means to calculate coefficient of variation
                                        "cell" means to calculate probability of each pixel being inundated
                filter_rate_func:
                (float or int)
                                        The rate at which the depth values will be ignored
    -----------

    -----------
    Returns:
                new_dataframe:
                (pandas dataframe)
                                        A dataframe includes x, y coordinates and filtered water depth values
    -----------

    """
    # Remove geometry data (x, y coordinates)
    nogeo_data_func = dataset_func.drop(['x_coord', 'y_coord'], axis=1)

    # Calculate mean data
    mean_data_func = nogeo_data_func.copy()
    mean_data_func['mean'] = mean_data_func.mean(axis=1)

    # Copy mean and filter data
    copy_mean_data_func = mean_data_func.copy()

    # Removing pixels having lower than 'filter_rate_func' m water depth values
    copy_mean_data_func.loc[copy_mean_data_func['mean'] < filter_rate_func, ['mean']] = -999

    # Manipulate and calculate necessary information
    if calculation_option == "mean":
        filter_data_func = copy_mean_data_func['mean']

    elif calculation_option == "sd":
        # Calculate standard deviation
        sd_data_func = nogeo_data_func.copy()
        sd_data_func['sd'] = sd_data_func.std(axis=1)

        # Copy standard deviation and filter data
        copy_sd_data_func = sd_data_func.copy()
        copy_sd_data_func.loc[copy_mean_data_func['mean'] == -999, ['sd']] = -999
        filter_data_func = copy_sd_data_func['sd']

    elif calculation_option == "cv":
        # Calculate coefficient of variation
        cv_data_func = nogeo_data_func.copy()
        cv_data_func['mean'] = cv_data_func.mean(axis=1)
        cv_data_func['sd'] = cv_data_func.std(axis=1)
        cv_data_func['cv'] = cv_data_func['sd'] / cv_data_func['mean'] * 100

        # Copy coefficient of variation and filter data
        copy_cv_data_func = cv_data_func.copy()
        copy_cv_data_func.loc[copy_mean_data_func['mean'] == -999, ['cv']] = -999
        filter_data_func = copy_cv_data_func['cv']

    else:
        # Calculate probability of each location getting inundated
        cell_data_func = nogeo_data_func.copy()
        cell_data_func['cell'] = (nogeo_data_func.shape[1] - (nogeo_data_func <= filter_rate_func).sum(axis=1)) / \
                                  nogeo_data_func.shape[1] * 100

        # Filter data
        cell_data_func.loc[copy_mean_data_func['mean'] == -999, ['cell']] = -999
        filter_data_func = cell_data_func['cell']

    # Create new dataframe with calculated information
    new_database = {'x': dataset_func['x_coord'],
                    'y': dataset_func['y_coord'],
                    f"{calculation_option}": filter_data_func}
    new_dataframe = pd.DataFrame(new_database)

    # Return new dataframe
    return new_dataframe


def area_calculation(dataset_func, resolution_func, filter_rate_func):
    """This function is to calculate area of each simulation

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                dataset_func:
                (pandas dataframe)
                                        A full dataset including necessary information
                resolution_func:
                (int or float)
                                        Resolution value in meter
                filter_rate_func:
                (float or int)
                                        The rate at which the depth values will be ignored
    -----------

    -----------
    Returns:
                area_dataframe:
                (pandas dataframe)
                                        Dataframe of simulations' areas
    -----------
    
    """
    # Remove geometry data (x, y coordinates)
    nogeo_data_func = dataset_func.drop(['x_coord', 'y_coord'], axis=1)

    # Get area dictionary
    area_dict = {}

    # Get area of each simulation
    for each_col in range(len(nogeo_data_func.columns)):
        # Get column num of each simulation
        column = nogeo_data_func[nogeo_data_func.columns[each_col]]

        # Get quantity of flooded cells
        num_flooded_cell = len(column[column > filter_rate_func])

        # Calculate area of each simulation
        area_func = num_flooded_cell * resolution_func ** 2

        # Add each area to area dictionary
        area_dict[f"{nogeo_data_func.columns[each_col]}"] = area_func

    # Write dictionary into dataframe
    area_dataframe = pd.DataFrame(data=area_dict, index=[0])

    return area_dataframe