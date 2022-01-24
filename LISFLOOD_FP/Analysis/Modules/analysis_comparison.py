# Prepare packages -----------------------------------------------------------------------------------------------------

import seaborn as sns                           # For plotting
import pandas as pd                             # For reading data
import numpy as np                              # For handling data (using clip function)

from runStatistic import calculation_dict       # For generating statistical dictionary

# ----------------------------------------------------------------------------------------------------------------------

def get_datalist(list_filename, transformation_selection, resolution_func_list, filter_rate_func):
    """This function is to call data from csv, filter by a specific value and store in pandas dataframe

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                list_filename:
                (list)
                                            A list of directory of data
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                resolution_func_list:
                (list)
                                            A list of resolution values in meters
                filter_rate_func:
                (float or int)
                                            The rate at which the depth values will be ignored
    -----------

    -----------
    Returns:
                data_list:
                (list)
                                            A list contains dictionaries of data
    -----------

    """
    # Create a list contains data
    data_list = []

    # Loop to get dataframe of each type of data
    for i in range(len(list_filename)):
        # Get data from csv file
        data_df = pd.read_csv(list_filename[i])

        if len(resolution_func_list) != 1:
            # Get data if there are many resolution values
            stat_data = calculation_dict(
                transformation_selection, data_df,
                resolution_func_list[i],
                filter_rate_func
            )
        else:
            # Get data if ther is only one resolution
            stat_data = calculation_dict(
                transformation_selection, data_df,
                resolution_func_list[0],
                filter_rate_func
            )

        # Add into list
        data_list.append(stat_data)

    return data_list


def statistic_dataframe(list_dataname, list_dataframe, calculation_option, filter_value):
    """This function is to create dataframe for statistical selection
    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                list_dataname:
                (list)
                                        A list contains name of data
                list_dataframe:
                (list)
                                        A list contains dataframe according to name of data
                calculation_option:
                (string)
                                        "mean" means to calculate mean
                                        "sd" means to calculate standard deviation
                                        "cv" means to calculate coefficient of variation
                                        "cell" means to calculate probability of each pixel being inundated
                filter_value:
                (int or float)
                                        Value to be filtered
    -----------

    -----------
    Returns:
                comparison_dictionary:
                (dictionary)
                                        A dictionary contains all dataframes of type of data
    -----------

    """
    # Get empty dictionary
    comparison_dictionary = dict()

    # Add keys and values into dictionary
    for i in range(len(list_dataname)):
        if calculation_option != "area":
            comparison_dictionary[f'{list_dataname[i]}'] = list_dataframe[i][calculation_option][
                calculation_option][list_dataframe[i][calculation_option][calculation_option] != filter_value]
        else:
            comparison_dictionary[f'{list_dataname[i]}'] = list_dataframe[i][calculation_option].iloc[0]

    return comparison_dictionary


def statistic_df_dict(list_dataname, list_dataframe, filter_value):
    """This function is to create a dictionary contains all dataframes of statistical selections

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                list_dataname:
                (list)
                                            A list contains name of data
                list_dataframe:
                (list)
                                            A list contains dataframe according to name of data
                filter_value:
                (int or float)
                                            Value to be filtered
    -----------

    -----------
    Returns:
                statistic_df_dictionary:
                (dictionary)
                                            A dictionary contains all dataframes of types of data
    -----------

    """
    # Get an empty dictionary
    statistic_df_dictionary = dict()

    # Statistical selections
    stat_selection = ['mean', 'sd', 'cv', 'cell', 'area']

    # Generate statistical selection and add into a dictionary
    for each_option in stat_selection:
        stat_selection_df = pd.DataFrame(statistic_dataframe(
            list_dataname,
            list_dataframe,
            each_option,
            filter_value
        ))

        statistic_df_dictionary[f'{each_option}'] = stat_selection_df

    return statistic_df_dictionary


def boxplots(statistic_df_dictionary, axis_func, title_func, y_label_func, calculation_option):
    """This function is to generate box plots

    -----------
    References:
                https://cduvallet.github.io/posts/2018/03/boxplots-in-python
                https://seaborn.pydata.org/generated/seaborn.boxplot.html
                https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.boxplot.html#matplotlib.axes.Axes.boxplot
                https://matplotlib.org/stable/gallery/pyplots/boxplot_demo_pyplot.html#sphx-glr-gallery-pyplots-boxplot-demo-pyplot-py
                https://stackoverflow.com/questions/35131798/tweaking-seaborn-boxplot
    -----------

    -----------
    Arguments:
                statistic_df_dictionary:
                (dictionary)
                                            A dictionary contains all dataframes of types of data
                axis_func:
                (axis in matplotlib)
                                            The ordinal number of plot in subplot
                title_func:
                (string)
                                            Title of side-by-side box plots
                y_label_func:
                (string)
                                            Label of y axis
                calculation_option:
                (string)
                                            "mean" means to calculate mean
                                            "sd" means to calculate standard deviation
                                            "cv" means to calculate coefficient of variation
                                            "cell" means to calculate probability of each pixel being inundated
    -----------

    -----------
    Returns:
                None.
    -----------

    """

    # Design mean visualisation
    meanprops = dict(
        marker='o',
        markerfacecolor='fuchsia',
        markersize=7,
        markeredgecolor='fuchsia'
    )

    # Design median visualisation
    medianprops = dict(
        color='aqua',
        linewidth=3
    )

    # Design outliers
    flierprops = dict(
        marker='o',
        markerfacecolor='0.75',
        markersize=5,
        linestyle='none'
    )

    # Design title
    title = f"Distributions of {title_func}"

    # Design x labels
    if calculation_option == 'mean':
        x_label = "Mean (m)"
    elif calculation_option == 'sd':
        x_label = "Standard deviation (m)"
    elif calculation_option == 'cv':
        x_label = "Coefficient of variation (%)"
    elif calculation_option == 'cell':
        x_label = "Proportion (%)"
    else:
        x_label = r"Areas ($m^{2}$)"

    # Draw boxplots
    sns.boxplot(data=statistic_df_dictionary[calculation_option],
                orient="h", palette='dark',
                showmeans=True, meanprops=meanprops,
                flierprops=flierprops,
                medianprops=medianprops,
                ax=axis_func)

    # Generate title and labels
    axis_func.set_title(title, pad=25, fontsize=25, fontweight='bold')
    axis_func.set_xlabel(x_label, fontsize=20, labelpad=38)
    axis_func.set_ylabel(y_label_func, rotation=-270, fontsize=20, labelpad=38)

    # Design size and style for ticks and labels
    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)
    axis_func.tick_params(direction='out', length=8, pad=10)

    # Design for 1e6
    axis_func.xaxis.offsetText.set_fontsize(12)



def kdeplots(statistic_df_dictionary,
             axis_func, title_func, x_range_func,
             legend_label_func,
             position_func, calculation_option,
             color_list_func, x_limit=""):
    """This function is to generate kernel density estimate (KDE). KDE represents the data using a continuous
    probability density curve in one or more dimensions. For more information, please visit
    https://seaborn.pydata.org/generated/seaborn.kdeplot.html

    -----------
    References:
                https://seaborn.pydata.org/generated/seaborn.kdeplot.html
                https://stackoverflow.com/questions/26218704/matplotlib-histogram-with-collection-bin-for-high-values
                https://stackoverflow.com/questions/37079573/does-seaborn-distplot-not-support-a-range
                https://stackoverflow.com/questions/45911709/limit-the-range-of-x-in-seaborn-distplot-kde-estimation

                https://stackoverflow.com/questions/44880444/how-to-increase-the-font-size-of-the-legend-in-my-seaborn-plot
                https://stackoverflow.com/questions/12402561/how-to-set-font-size-of-matplotlib-axis-legend

                https://stats.stackexchange.com/questions/362798/can-pmf-have-value-greater-than-1
                https://github.com/mwaskom/seaborn/issues/479
                
    -----------

    -----------
    Arguments:
                statistic_df_dictionary:
                (dictionary)
                                            A dictionary contains all dataframes of types of data
                axis_func:
                (axis in matplotlib)
                                            The ordinal number of plot in subplot
                title_func:
                (string)
                                            Title of side-by-side box plots
                x_range_func:
                (list)
                                            A list contains min, max, step in order
                legend_label_func:
                (string)
                                            Title of legends
                position_func:
                (string)
                                            Name of position to locate legend. Please vistit here for more information
                                            
                calculation_option:
                (string)
                                            "mean" means to calculate mean
                                            "sd" means to calculate standard deviation
                                            "cv" means to calculate coefficient of variation
                                            "cell" means to calculate probability of each pixel being inundated
                color_list_func:
                (list)
                                            A list contains color names
                x_limit:
                (string)
                                            Additional information for x tick labels
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    
    # Design title
    title = f"Distributions of {title_func}"

    # Design x labels
    if calculation_option == 'mean':
        x_label = "Mean (m)"
    elif calculation_option == 'sd':
        x_label = "Standard deviation (m)"
    elif calculation_option == 'cv':
        x_label = "Coefficient of variation (%)"
    elif calculation_option == 'cell':
        x_label = "Proportion (%)"
    else:
        x_label = r"Areas ($m^{2}$)"

    # Get column names from selected dataframes
    col_names = statistic_df_dictionary[calculation_option].columns

    # Draw kde plot
    for num_name in range(len(color_list_func)):
        sns.kdeplot(
            np.clip(
                statistic_df_dictionary[calculation_option][statistic_df_dictionary[calculation_option].columns[
                    num_name]],
                x_range_func[0], x_range_func[1]
            ),
            fill=True, linewidth=2,
            label=col_names[num_name],
            clip=(x_range_func[0], x_range_func[1]),
            alpha=0.2,
            color=color_list_func[num_name],
            ax=axis_func
        )

    # Add legends
    axis_func.legend(prop={'size': 20}, title=legend_label_func, title_fontsize=18, loc=position_func)

    # Generate title and labels
    axis_func.set_title(title, pad=25, fontsize=25, fontweight='bold')
    axis_func.set_xlabel(x_label, fontsize=20, labelpad=38)
    axis_func.set_ylabel("Density", rotation=-270, fontsize=20, labelpad=38)

    # Remove grid background lines (including x, y lines)
    axis_func.grid(False)
    axis_func.spines['top'].set_visible(False)
    axis_func.spines['right'].set_visible(False)
    axis_func.spines['bottom'].set_visible(False)
    axis_func.spines['left'].set_visible(False)

    ## Set x label, for more information, please visit
    # https://stackoverflow.com/questions/26218704/matplotlib-histogram-with-collection-bin-for-high-values
    # For dtype: https://numpy.org/doc/stable/reference/arrays.dtypes.html

    # Set tick
    x_range = np.arange(x_range_func[0], x_range_func[1], x_range_func[2])
    xlabel_arr = np.array(np.round(x_range[:], 2), dtype='float')
    axis_func.set_xticks(xlabel_arr)

    # Set tick label
    xlabel_arr_cop = np.array(np.round(x_range[:], 2), dtype='str')
    xlabel_arr_cop[-1] = '{0:.2f}{1}'.format(xlabel_arr[-1], x_limit)
    axis_func.set_xticklabels(xlabel_arr_cop)

    # Design size and style for ticks and labels
    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)
    axis_func.tick_params(direction='out', length=8, pad=10)
    

