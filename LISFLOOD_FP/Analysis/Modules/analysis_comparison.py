# Prepare packages -----------------------------------------------------------------------------------------------------

import seaborn as sns                           # For plotting
import pandas as pd                             # For reading data
import numpy as np                              # For handling data (using clip function)
import matplotlib


from runStatistic import calculation_dict       # For generating statistical dictionary

# ----------------------------------------------------------------------------------------------------------------------

def get_datalist(list_filename, transformation_selection, resolution_func_list, filter_rate_func, building_path,
                 dif_path_list, raster=True,
                 rectangle=True, switch=False):
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
                building_path:
                (string)
                                            Path of file containing building polygons
                dif_path_list:
                (list)
                                            A list of paths of different file versions
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
                filter_rate_func,
                building_path,
                dif_path_list[i],
                raster,
                rectangle,
                switch
            )
        else:
            # Get data if there is only one resolution
            stat_data = calculation_dict(
                transformation_selection, data_df,
                resolution_func_list[0],
                filter_rate_func,
                building_path,
                dif_path_list[i],
                raster,
                rectangle,
                switch
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
        if calculation_option not in ["area", "building"]:
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
    stat_selection = ['mean', 'sd', 'cv', 'cell', 'area', 'building']

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


def boxplots(statistic_df_dictionary, axis_func, title_func, y_label_func, calculation_option, color_dict_func,
             reverse_order=False):
    """This function is to generate box plots

    -----------
    References:
                https://cduvallet.github.io/posts/2018/03/boxplots-in-python
                https://seaborn.pydata.org/generated/seaborn.boxplot.html
                https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.boxplot.html#matplotlib.axes.Axes.boxplot
                https://matplotlib.org/stable/gallery/pyplots/boxplot_demo_pyplot.html#sphx-glr-gallery-pyplots-boxplot-demo-pyplot-py
                https://stackoverflow.com/questions/35131798/tweaking-seaborn-boxplot

                https://www.python-graph-gallery.com/33-control-colors-of-boxplot-seaborn
                https://www.python-graph-gallery.com/32-custom-boxplot-appearance-seaborn
                https://stackoverflow.com/questions/46735745/how-to-control-scientific-notation-in-matplotlib
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
                color_dict_func:
                (dictionary)
                                            A dictionary contains color names
    -----------

    -----------
    Returns:
                None.
    -----------

    """

    # Design mean visualisation
    meanprops = dict(
        marker='o',
        markerfacecolor='lime',
        markersize=7,
        markeredgecolor='darkgreen'
    )

    # Design median visualisation
    medianprops = dict(
        color='aqua',
        linewidth=1.8
    )

    # Design outliers
    flierprops = dict(
        marker='o',
        markersize=8,
        markerfacecolor='none',
        markeredgecolor='red',
        linewidth=0.2,
        linestyle='none'
    )

    # # Design whisker
    # whiskerprops = dict(color='black', linewidth=1)

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
    elif calculation_option == 'area':
        x_label = r'Areas (x100 $\mathrm{m}^2$)'
    else:
        x_label = "Number of buildings"

    # Check if reverse_order is in use
    if reverse_order:
        dat = statistic_df_dictionary[calculation_option]
        dat_reversed = dat[dat.columns[::-1]]
    else:
        dat_reversed = statistic_df_dictionary[calculation_option]

    # Draw boxplots
    sns.boxplot(data=dat_reversed,
                orient="h",
                showmeans=True,
                # whiskerprops=whiskerprops,
                linewidth=1,
                meanprops=meanprops,
                flierprops=flierprops,
                medianprops=medianprops,
                palette=color_dict_func,
                ax=axis_func)

    # Generate title and labels
    axis_func.set_title(title, pad=25, fontsize=25, fontweight='bold')
    axis_func.set_xlabel(x_label, fontsize=20, labelpad=38)
    axis_func.set_ylabel(y_label_func, rotation=-270, fontsize=20, labelpad=38)

    # Control scientific notation
    # Refer here for more information
    # https://stackoverflow.com/questions/46735745/how-to-control-scientific-notation-in-matplotlib
    if calculation_option == 'area':
        axis_func.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x/100))))
    else:
        pass

    # Design size and style for ticks and labels
    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)
    axis_func.tick_params(direction='out', length=8, pad=10)

def skewness(series):
    return (((series - series.mean()) / series.std(ddof=0)) ** 3).mean()

def kdeplots(statistic_df_dictionary,
             axis_func, title_func, x_range_func,
             legend_label_func,
             position_func, calculation_option,
             color_list_func, x_limit="", reverse_legend=False):
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
                                            Name of position to locate legend. Please visit here for more information
                                            
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
                reverse_legend:
                (boolean)
                                            To reverse the order of variables
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    
    # Design title
    title = f"Distributions of {title_func}"

    # Get column names from selected dataframes
    col_names = statistic_df_dictionary[calculation_option].columns

    # Design labels
    if calculation_option == 'mean':
        x_label = "Mean (m)"
    elif calculation_option == 'sd':
        x_label = "Standard deviation (m)"
    elif calculation_option == 'cv':
        x_label = "Coefficient of variation (%)"
    else:
        x_label = "Proportion (%)"

    # Draw kde plot
    for num_name in range(len(color_list_func)):
        # Check to see if reverse_legend is in use
        if reverse_legend:
            label_legend = col_names[::-1]
            color_each = color_list_func[::-1]
            dat = statistic_df_dictionary[calculation_option]
            dat_reversed = dat[dat.columns[::-1]]
        else:
            label_legend = col_names
            color_each = color_list_func
            dat_reversed = statistic_df_dictionary[calculation_option]

        # Draw plot
        sns.kdeplot(
            np.clip(
                dat_reversed[dat_reversed.columns[num_name]],
                x_range_func[0], x_range_func[1]
            ),
            fill=True, linewidth=1.8,
            label=label_legend[num_name],
            clip=(x_range_func[0], x_range_func[1]),
            alpha=0.25,
            color=color_each[num_name],
            ax=axis_func
        )

        print(skewness(dat_reversed[dat_reversed.columns[num_name]]), "-", label_legend[num_name])

    # Add legends
    axis_func.legend(prop={'size': 14}, title=legend_label_func,
                     title_fontproperties={'size': 17, 'weight': 'bold'},
                     frameon=False,
                     loc=position_func)

    # Generate title and labels
    axis_func.set_title(title, pad=25, fontsize=25, fontweight='bold')
    axis_func.set_xlabel(x_label, fontsize=20, labelpad=38)
    axis_func.set_ylabel("Probability density", rotation=-270, fontsize=20, labelpad=38)

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
    # Refer here for information about signs of larger/less or equal to
    # https://matplotlib.org/stable/tutorials/text/mathtext.html
    # https://matplotlib.org/2.0.2/users/mathtext.html
    # https://stackoverflow.com/questions/18428823/python-matplotlib-less-than-or-equal-to-symbol-in-text
    # https://stackoverflow.com/questions/48646669/make-matplotlib-display-serif-font-in-math-mode
    if calculation_option != 'cell':
        xlabel_arr_cop = np.array(np.round(x_range[:], 2), dtype='str')
        xlabel_arr_cop[-1] = '{1}{0:.2f}'.format(xlabel_arr[-1], x_limit)
        axis_func.set_xticklabels(xlabel_arr_cop)
    else:
        xlabel_arr_cop = np.array(np.round(x_range[:], 2), dtype='str')
        xlabel_arr_cop[0] = '{1}{0:.2f}'.format(xlabel_arr[0], x_limit)
        axis_func.set_xticklabels(xlabel_arr_cop)

    # Design size and style for ticks and labels
    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)
    axis_func.tick_params(direction='out', length=8, pad=10)


def kdeplots_vers2(statistic_df_dictionary,
                   axis_func, title_func,
                   legend_label_func,
                   position_func, calculation_option,
                   color_list_func, reverse_legend=False):
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
                legend_label_func:
                (string)
                                            Title of legends
                position_func:
                (string)
                                            Name of position to locate legend. Please visit here for more information

                calculation_option:
                (string)
                                            "mean" means to calculate mean
                                            "sd" means to calculate standard deviation
                                            "cv" means to calculate coefficient of variation
                                            "cell" means to calculate probability of each pixel being inundated
                color_list_func:
                (list)
                                            A list contains color names
                reverse_legend:
                (boolean)
                                            To reverse the order of variables
    -----------

    -----------
    Returns:
                None.
    -----------

    """

    # Design title
    title = f"Distributions of {title_func}"

    # Get column names from selected dataframes
    col_names = statistic_df_dictionary[calculation_option].columns

    # Design labels
    if calculation_option == 'area':
        x_label = f'Areas (x100 '
        x_label += r'$\mathrm{m}^2$)'
    else:
        x_label = "Number of buildings"

    # Draw kde plot
    for num_name in range(len(color_list_func)):
        # Check to see if reverse_legend is in use
        if reverse_legend:
            label_legend = col_names[::-1]
            color_each = color_list_func[::-1]
            dat = statistic_df_dictionary[calculation_option]
            dat_reversed = dat[dat.columns[::-1]]
        else:
            label_legend = col_names
            color_each = color_list_func
            dat_reversed = statistic_df_dictionary[calculation_option]

        sns.kdeplot(
            dat_reversed[dat_reversed.columns[num_name]],
            fill=True, linewidth=1.8,
            label=label_legend[num_name],
            alpha=0.3,
            color=color_each[num_name],
            ax=axis_func
        )

    # Add legends
    # Refer here for more information
    # https://stackoverflow.com/questions/25540259/remove-or-adapt-border-of-frame-of-legend-using-matplotlib
    # https://stackoverflow.com/questions/33661080/python-matplotlib-plot-text-will-not-align-left
    # https://matplotlib.org/stable/gallery/text_labels_and_annotations/text_alignment.html
    axis_func.legend(prop={'size': 14}, title=legend_label_func,
                     title_fontproperties={'size': 17, 'weight': 'bold'},
                     frameon=False,
                     loc=position_func)

    # Generate title and labels
    axis_func.set_title(title, pad=25, fontsize=25, fontweight='bold')
    axis_func.set_xlabel(x_label, fontsize=20, labelpad=38)
    axis_func.set_ylabel("Probability density", rotation=-270, fontsize=20, labelpad=38)

    # Remove grid background lines (including x, y lines)
    axis_func.grid(False)
    axis_func.spines['top'].set_visible(False)
    axis_func.spines['right'].set_visible(False)
    axis_func.spines['bottom'].set_visible(False)
    axis_func.spines['left'].set_visible(False)

    # Remove scientific notation
    if calculation_option == 'area':
        axis_func.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x / 100))))
    else:
        pass
    
    # Design size and style for ticks and labels
    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)
    axis_func.tick_params(direction='out', length=8, pad=10)



    
    
    

