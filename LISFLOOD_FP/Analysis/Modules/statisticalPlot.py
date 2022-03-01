# Prepare packages -----------------------------------------------------------------------------------------------------
import matplotlib.ticker
import numpy as np                                              # For data manipulation
import seaborn as sns                                           # For plotting density and frequency
from mpl_toolkits.axes_grid1 import make_axes_locatable         # For vertical color bar
from colorDevelop import get_gradient_cmap                      # For gradient color development

# ----------------------------------------------------------------------------------------------------------------------


# STATISTICAL PLOT #####################################################################################################
def round_down(number_func, decimals=0):
    """This function is to round down the number.
    This code is based on the link in reference (by Priyankur Sarkar)

    -----------
    References: https://www.knowledgehut.com/blog/programming/python-rounding-numbers
    -----------

    -----------
    Arguments:
                number_func:
                (float or int)
                                        Number needs rounding down
                decimals:
                (int)
                                        A specific decimal that the number can be rounded down to
    -----------

    -----------
    Returns:
                (float or int)
                                        Number that is rounded down
    -----------

    """
    return np.floor(number_func * (10 ** decimals)) / (10 ** decimals)


def plotting_map(filtered_data_func,
                 resolution_func, calculation_option,
                 axis_func, fig_func,
                 hex_list_func, colorbar_position,
                 extend_colorbar=None,
                 color_range=None,
                 add_title="",):
    """This function is to plot water depth on basemap

    -----------
    References: https://stackoverflow.com/questions/41897544/make-a-contour-plot-by-using-three-1d-arrays-in-python
                https://stackoverflow.com/questions/44669616/contour-in-matplotlib-does-not-plot-specified-number-of-contours
                https://stackoverflow.com/questions/48487346/filled-contour-using-class-labels
                https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/tricontour_smooth_user.html#sphx-glr-gallery-images-contours-and-fields-tricontour-smooth-user-py
                https://www.python-course.eu/matplotlib_contour_plot.php
                https://jakevdp.github.io/PythonDataScienceHandbook/04.04-density-and-contour-plots.html
                https://www.earthdatascience.org/tutorials/visualize-digital-elevation-model-contours-matplotlib/
                https://alex.miller.im/posts/contour-plots-in-python-matplotlib-x-y-z/

                https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.contour.html
                https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.tricontourf.html
                https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.tricontour.html#matplotlib.pyplot.tricontour

                https://contextily.readthedocs.io/en/latest/providers_deepdive.html
                https://stackoverflow.com/questions/3777861/setting-y-axis-limit-in-matplotlib
                https://www.tutorialspoint.com/how-do-i-adjust-offset-the-colorbar-title-in-matplotlib
                https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph

                https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
                https://kaleidoscopicdiaries.wordpress.com/2015/05/30/distance-between-axes-label-and-axes-in-matplotlib/
                https://geopandas.org/en/latest/gallery/matplotlib_scalebar.html
                https://stackoverflow.com/questions/8263769/hide-contour-linestroke-on-pyplot-contourf-to-get-only-fills

                https://matplotlib.org/stable/gallery/axes_grid1/demo_colorbar_with_axes_divider.html
                https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    -----------

    -----------
    Arguments:
                filtered_data_func:
                (pandas dataframe)
                                            Filtered dataset contains nodata values as -999
                resolution_func:
                (int or float)
                                            Resolution value in meter
                calculation_option:
                (string)
                                            "mean" means to calculate mean
                                            "sd" means to calculate standard deviation
                                            "cv" means to calculate coefficient of variation
                                            "cell" means to calculate probability of each pixel being inundated
                fig_func:
                (figure in matplotlib)
                                            Figure from matplotlib subplot
                axis_func:
                (axis in matplotlib)
                                            Subplot ordinal number of the contour map on the big plot
                hex_list_func:
                                            A list of hex code colors (including or not including #)
                colorbar_position:
                                            "horizontal" or "vertical"
                extend_colorbar:
                                            Whether having or not having arrow to indicate
                                            the larger or smaller values in colorbar.
                                            Selections are "neither", "both", "min", "max".
                                            For more information, please visit
                                            https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.tricontourf.html
                color_range:
                (list)
                                            A list includes min and max to be displayed on the color scale
                add_title:
                (string)
                                            Add more words into title
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Titles and labels:
    if calculation_option == 'cell':
        # Title for colorbar
        name_map = "Proportion of simulations\nof each area being inundated,"
        name_map += f"\nresolution = {resolution_func} meters"
        name_map += add_title

        # Title for contour map
        axis_func.set_title(name_map, pad=25, fontsize=25, fontweight='bold')
        axis_func.set_xlabel("NZTM, east (m)", fontsize=20, labelpad=38)
        axis_func.set_ylabel("NZTM, north (m)", fontsize=20, labelpad=38, rotation=-270)

        # Max level
        if color_range is None:
            min_map_level = 0
            max_map_level = 102
            map_level = np.arange(min_map_level, max_map_level, 1)
        else:
            min_map_level = color_range[0]
            max_map_level = color_range[1]
            map_level = np.arange(min_map_level, max_map_level, color_range[2])

    else:
        # Title for colorbar
        if calculation_option == 'sd':
            name_map = f"Standard deviation of water depth,\nresolution = {resolution_func} meters"
            name_map += add_title
        elif calculation_option == 'cv':
            name_map = f"Coefficient of variation of water depth,\nresolution = {resolution_func} meters"
            name_map += add_title
        else:
            name_map = f"Mean of water depth,\nresolution = {resolution_func} meters"
            name_map += add_title

        # Title for contour map
        axis_func.set_title(f"{name_map}", pad=25, fontsize=25, fontweight='bold')
        axis_func.set_xlabel("NZTM, east (m)", fontsize=20, labelpad=38)
        axis_func.set_ylabel("NZTM, north (m)", fontsize=20, labelpad=38, rotation=-270)

        # Level information
        if color_range is None:
            min_map_level = round_down(
                np.min(filtered_data_func[f'{calculation_option}'][filtered_data_func[f'{calculation_option}'] != -999]), 3)
            max_map_level = np.max(
                filtered_data_func[f'{calculation_option}'][filtered_data_func[f'{calculation_option}'] != -999]) + 0.000001
            map_level = np.arange(min_map_level, max_map_level, 0.1)
        else:
            min_map_level = color_range[0]
            max_map_level = color_range[1]
            map_level = np.arange(min_map_level, max_map_level, color_range[2])

    # x, y, z coordinates information
    x_func = filtered_data_func['x']
    y_func = filtered_data_func['y']
    z_func = filtered_data_func[f'{calculation_option}']

    # Build up contour map
    contour_map_func = axis_func.tricontourf(x_func, y_func, z_func, levels=map_level,
                                             cmap=get_gradient_cmap(hex_list_func),
                                             alpha=1, antialiased=True, extend=extend_colorbar)

    # Improve visualisation of map (x, y) axes
    axis_func.tick_params(direction='out', length=8, pad=10)  # For x, y axes' ticks
    axis_func.yaxis.offsetText.set_fontsize(
        12)  # For 1e6 if there is no axis.ticklabel_format(useOffset=False, style='plain')
    axis_func.xaxis.offsetText.set_fontsize(
        12)  # For 1e6 if there is no axis.ticklabel_format(useOffset=False, style='plain')

    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)

    # Remove grid background lines (including x, y lines)
    axis_func.grid(False)
    axis_func.spines['top'].set_visible(False)
    axis_func.spines['right'].set_visible(False)
    axis_func.spines['bottom'].set_visible(False)
    axis_func.spines['left'].set_visible(False)

    # Remove whiteness
    for color_line in contour_map_func.collections:
        color_line.set_edgecolor("face")

    #     # Call the map
    #     ctx.add_basemap(ax=axis_func, crs=2193,
    #                     source=ctx.providers.OpenTopoMap,
    #                     attribution_size=10, zoom=17)

    #     # Scale bar
    #     map_scalebar = ScaleBar(10, font_properties={'weight': 'bold', 'size': 15},
    #                             pad=0.5,
    #                             box_color=None,
    #                             box_alpha=0,
    #                             color='k',
    #                             scale_formatter=lambda value, unit: f'{value} {unit}')

    #     axis_func.add_artist(map_scalebar)

    # Colorbar
    if colorbar_position == 'horizontal':
        # Set up colorbar position
        divider = make_axes_locatable(axis_func)
        cax = divider.append_axes("bottom", size="3%", pad=1.8)

        map_colorbar = fig_func.colorbar(contour_map_func, orientation=colorbar_position, cax=cax, ax=axis_func)
        map_colorbar.ax.tick_params(axis='x', direction='out', length=6, pad=10, labelsize=15)

        # Set up colorbar title
        if calculation_option == "cell":
            map_colorbar.set_label(f"Proportion (%)", labelpad=30, fontsize=20)
        elif calculation_option == "cv":
            map_colorbar.set_label(f"Coefficient of variation (%)", labelpad=30, fontsize=20)
        elif calculation_option == "sd":
            map_colorbar.set_label(f"Standard deviation (m)", labelpad=30, fontsize=20)
        else:
            map_colorbar.set_label(f"{calculation_option.capitalize()} (m)", labelpad=30, fontsize=20)


    else:
        # Set up colorbar position
        cax = fig_func.add_axes([axis_func.get_position().x1 + 0.02,
                                axis_func.get_position().y0, 0.02,
                                axis_func.get_position().height])

        map_colorbar = fig_func.colorbar(contour_map_func, cax=cax, pad=5, ax=axis_func)
        map_colorbar.ax.tick_params(axis='y', direction='out', length=6, pad=10, labelsize=15)

        # Set up colorbar title
        if calculation_option == "cell":
            map_colorbar.set_label(f"Proportion (%)", labelpad=30, fontsize=20)
        elif calculation_option == "cv":
            map_colorbar.set_label(f"Coefficient of variation (%)", labelpad=30, fontsize=20)
        elif calculation_option == "sd":
            map_colorbar.set_label(f"Standard deviation (m)", labelpad=30, fontsize=20)
        else:
            map_colorbar.set_label(f"{calculation_option.capitalize()} (m)", labelpad=30, fontsize=20)



def plotting_histogram(filtered_data_func,
                       resolution_func, calculation_option,
                       plt_func, axis_func, hex_list_func,
                       color_range,
                       num_bin_func,
                       text_box_loc,
                       x_limit_sign="",
                       y_limit=None,
                       add_title=""):
    """This function is to plot histogram of information regarding water depth

    -----------
    References: https://stackoverflow.com/questions/23061657/plot-histogram-with-colors-taken-from-colormap (first)
                https://stdworkflow.com/67/attributeerror-rectangle-object-has-no-property-normed-solution
                https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib
                https://stackoverflow.com/questions/63723514/userwarning-fixedformatter-should-only-be-used-together-with-fixedlocator
                https://stackoverflow.com/questions/38667728/matplotlib-histogram-scale-y-axis-by-a-constant-factor
                https://stackoverflow.com/questions/7559242/matplotlib-strings-as-labels-on-x-axis

                https://www.kite.com/python/answers/how-to-change-font-of-a-plot-using-matplotlib-in-python
    -----------

    -----------
    Arguments:
                filtered_data_func:
                (pandas dataframe)
                                        Filtered dataset contains nodata values as -999
                resolution_func:
                (int or float)
                                        Resolution value in meter
                calculation_option:
                (string)
                                        "mean" means to calculate mean
                                        "sd" means to calculate standard deviation
                                        "cv" means to calculate coefficient of variation
                                        "cell" means to calculate probability of each pixel being inundated
                plt_func:
                                        Short writing of matplotlib.pyplot
                axis_func:
                (axis in matplotlib)
                                        Subplot ordinal number of the contour map on the big plot
                hex_list_func:
                (string)
                                        A list of hex code colors (including or not including #)

                The two parameters below are for designing the range for histogram

                color_range:
                (list)
                                        A list includes min and max to be displayed on the color scale
                num_bin_func:
                (array or list)
                                        An array/ a list contains the range of bins
                text_box_loc:
                (list)
                                        A list of text box coordinates [x, y]
                x_limit_sign:
                (string)
                                        Additional information for x tick labels
                y_limit:
                (string)
                                        New limit of y range
                add_title:
                (string)
                                        Add more words into title
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Get z values
    z_func = filtered_data_func[f'{calculation_option}']

    # Assign data into a variable
    data_func = z_func[z_func != -999]

    # Design the range of histogram
    # Get max values
    if color_range[1] is None:
        max_map_level = np.max(
            filtered_data_func[f'{calculation_option}'][filtered_data_func[f'{calculation_option}'] != -999]) + 0.000001
    else:
        max_map_level = color_range[1]
    # Develop the range for x axis
    axis_func.set_xticks(np.arange(color_range[0], max_map_level, color_range[2]))

    # Get values of n, bins, and patches. Please visit here for more information
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.hist.html#matplotlib.axes.Axes.hist
    # Besides, to find out how to collect all values for one bin representing only one value,
    # please visit https://stackoverflow.com/questions/26218704/matplotlib-histogram-with-collection-bin-for-high-values
    n, bins, patches = axis_func.hist(np.clip(data_func, num_bin_func[0], num_bin_func[-1]),
                                      bins=num_bin_func,
                                      density=False)


    # Get values of columns (scale values to interval [0,1])
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    col = bin_centers - min(bin_centers)
    col /= max(col)

    # Get gradient color
    gradient_color = get_gradient_cmap(hex_list_func)

    # Loop to assign color and values into each patch
    for c, p in zip(col, patches):
        plt_func.setp(p, 'facecolor', gradient_color(c))

    # Set x label, for more information, please visit
    # https://stackoverflow.com/questions/26218704/matplotlib-histogram-with-collection-bin-for-high-values
    # For dtype: https://numpy.org/doc/stable/reference/arrays.dtypes.html
    x_range = np.arange(color_range[0], color_range[1], color_range[2])
    xlabel_arr = np.array(np.round(x_range[:], 2), dtype='str')
    xlabel_arr[-1] = f'{x_limit_sign}{xlabel_arr[-1]}'
    axis_func.set_xticklabels(xlabel_arr)

    # Remove grid background lines (including x, y lines)
    axis_func.grid(False)
    axis_func.spines['top'].set_visible(False)
    axis_func.spines['right'].set_visible(False)
    axis_func.spines['bottom'].set_visible(False)
    axis_func.spines['left'].set_visible(False)

    # Set titles and labels
    if calculation_option == 'cell':
        # Title for contour map and x label
        name_hist = "Histogram of proportion of simulations\nof each area being inundated,"
        name_hist += f"\nresolution = {resolution_func} meters"
        name_hist += add_title
        axis_func.set_xlabel("Proportion (%)", fontsize=20, labelpad=38)


        # Set up new y tick axis
        # Refer here for more information
        # https://stackoverflow.com/questions/63723514/userwarning-fixedformatter-should-only-be-used-together-with-fixedlocator
        # https://stackoverflow.com/questions/38667728/matplotlib-histogram-scale-y-axis-by-a-constant-factor
        # https://stackoverflow.com/questions/7559242/matplotlib-strings-as-labels-on-x-axis

        axis_func.set_ylim(0, y_limit)                                              # Set up new range
        num_tick = axis_func.get_yticks().tolist()                                  # Convert new range into list
        num_tick_label = [str(int(val * (resolution_func**2) / 100)) for val in num_tick] # Convert cells into areas
        num_tick_label[-2] = f'{num_tick_label[-2]}+'                               # Add upper limit (by adding +)
        axis_func.set_yticks(num_tick[:-1])                                         # Set up new tick sign
        axis_func.set_yticklabels(num_tick_label[:-1])                              # Set up new tick label

    elif calculation_option == "cv":
        # Title for contour map and x label
        name_hist = "Histogram of coefficient of variation of water depth,"
        name_hist += f"\nresolution = {resolution_func} meters"
        name_hist += add_title
        axis_func.set_xlabel('Coefficient of variation (%)', fontsize=20, labelpad=38)

        # Set new y tick axis
        if y_limit == None:
            num_tick = axis_func.get_yticks().tolist()                                     # Convert new range into list
            num_tick_label = [str(int(val * (resolution_func ** 2) / 100)) for val in num_tick]  # Convert cells into areas
            axis_func.set_yticks(num_tick)                                                 # Set up new tick sign
            axis_func.set_yticklabels(num_tick_label)                                      # Set up new tick label
        else:
            axis_func.set_ylim(0, y_limit)  # Set up new range
            num_tick = axis_func.get_yticks().tolist()  # Convert new range into list
            num_tick_label = [str(int(val * (resolution_func ** 2) / 100)) for val in
                              num_tick]  # Convert cells into areas
            num_tick_label[-2] = f'{num_tick_label[-2]}+'  # Add upper limit (by adding +)
            axis_func.set_yticks(num_tick[:-1])  # Set up new tick sign
            axis_func.set_yticklabels(num_tick_label[:-1])

    elif calculation_option == "sd":
        # Title for contour map and x label
        name_hist = "Histogram of standard deviation of water depth,"
        name_hist += f"\nresolution = {resolution_func} meters"
        name_hist += add_title
        axis_func.set_xlabel('Standard deviation (m)', fontsize=20, labelpad=38)

        # Set new y tick axis
        if y_limit == None:
            num_tick = axis_func.get_yticks().tolist()                                     # Convert new range into list
            num_tick_label = [str(int(val * (resolution_func ** 2) / 100)) for val in num_tick]  # Convert cells into areas
            axis_func.set_yticks(num_tick)                                                 # Set up new tick sign
            axis_func.set_yticklabels(num_tick_label)                                     # Set up new tick label
        else:
            axis_func.set_ylim(0, y_limit)  # Set up new range
            num_tick = axis_func.get_yticks().tolist()  # Convert new range into list
            num_tick_label = [str(int(val * (resolution_func ** 2) / 100)) for val in
                              num_tick]  # Convert cells into areas
            num_tick_label[-2] = f'{num_tick_label[-2]}+'  # Add upper limit (by adding +)
            axis_func.set_yticks(num_tick[:-1])  # Set up new tick sign
            axis_func.set_yticklabels(num_tick_label[:-1])

    else:
        # Title for contour map
        name_hist = f"Histogram of {calculation_option} of water depth,"
        name_hist += f"\nresolution = {resolution_func} meters"
        name_hist += add_title
        axis_func.set_xlabel(f"{calculation_option.capitalize()} (m)", fontsize=20, labelpad=38)
        
        # Set new y tick axis
        if y_limit == None:
            num_tick = axis_func.get_yticks().tolist()                                     # Convert new range into list
            num_tick_label = [str(int(val * (resolution_func ** 2) / 100)) for val in num_tick]  # Convert cells into areas
            axis_func.set_yticks(num_tick[:-1])                                            # Set up new tick sign
            axis_func.set_yticklabels(num_tick_label[:-1])                                # Set up new tick label
        else:
            axis_func.set_ylim(0, y_limit)  # Set up new range
            num_tick = axis_func.get_yticks().tolist()  # Convert new range into list
            num_tick_label = [str(int(val * (resolution_func ** 2) / 100)) for val in
                              num_tick]  # Convert cells into areas
            num_tick_label[-2] = f'{num_tick_label[-2]}+'  # Add upper limit (by adding +)
            axis_func.set_yticks(num_tick[:-1])  # Set up new tick sign
            axis_func.set_yticklabels(num_tick_label[:-1])

    # Ticks, title, and y label
    axis_func.tick_params(direction='out', length=8, pad=10)
    axis_func.set_title(name_hist, pad=25, fontsize=25, fontweight='bold')
    axis_func.set_ylabel(r'Areas (x100 $\mathrm{m}^2$)', fontsize=20, labelpad=38)

    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)

    # Create text string
    textstr = "mean = {0:.3f}".format(data_func.mean(axis=0))
    textstr += "\nstdev = {0:.3f}".format(data_func.std(axis=0))

    # place a text box in upper left in axes coords
    # Refer here for more information
    # https://stackoverflow.com/questions/50869424/increase-line-separation-in-matplotlib-annotation-text/50888491
    axis_func.text(text_box_loc[0], text_box_loc[1], textstr, transform=axis_func.transAxes, fontsize=13,
                   linespacing=1.5,
                   fontweight='normal',
                   fontstyle='italic',
                   horizontalalignment='left',
                   verticalalignment='top')

# END STATISTICAL PLOT #################################################################################################



# AREA PLOT ############################################################################################################
def plot_area(axis_func, area_dataframe_func, resolution_func, text_box_loc,
              density='y', add_title=""):
    """This function is to plot histogram of areas of simulations

    -----------
    References:
                https://stackoverflow.com/questions/69524514/how-to-modify-the-kernel-density-estimate-line-in-a-sns-histplot
                https://seaborn.pydata.org/generated/seaborn.distplot.html

                https://stackoverflow.com/questions/65400669/how-to-generate-two-separate-y-axes-for-a-histogram-on-the-same-figure-in-seabor
                https://stackoverflow.com/questions/26752464/how-do-i-align-gridlines-for-two-y-axis-scales-using-matplotlib
                https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.hist.html#matplotlib.axes.Axes.hist
                https://seaborn.pydata.org/generated/seaborn.kdeplot.html

                https://seaborn.pydata.org/generated/seaborn.histplot.html#seaborn.histplot
                https://seaborn.pydata.org/generated/seaborn.distplot.html
                https://stackoverflow.com/questions/27671748/how-to-print-y-axis-label-horizontally-in-a-matplotlib-pylab-chart
                https://stackoverflow.com/questions/24391892/printing-subscript-in-python

                https://stackoverflow.com/questions/45037386/trouble-aligning-ticks-for-matplotlib-twinx-axes (best
                answer for align two axis)
                https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib
    -----------

    -----------
    Arguments:
                axis_func:
                (axis in matplotlib)
                                        The subplot
                area_dataframe_func:
                (pandas dataframe)
                                        Dataframe of simulations' areas
                resolution_func:
                (int or float)
                                        Resolution value in meter
                text_box_loc:
                (list)
                                        A list of text box coordinates [x, y]
                add_title:
                (string)
                                        Add more words into title
                density:
                (string)
                                        Performing "density" of "frequency"
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Get areas' values under array format
    area_values = area_dataframe_func.iloc[0].to_numpy()

    # # Get unit of area
    # unit = np.power(resolution_func, 2)

    # Number of bins/simulations
    num_bin = len(area_values)

    # Get histogram and frequency
    if density == 'y':
        # Draw kde (kernel density estimate) plot
        sns.kdeplot(area_values / 100, color='crimson', ax=axis_func)

        # Set y label
        axis_func.set_ylabel("Probability density", fontsize=15, labelpad=38)


    elif density == 'n':
        # Draw histogram
        # sns.distplot(area_values/unit, bins=num_bin, kde=False, ax=axis_func,
        #              hist_kws={"histtype": "bar", "linewidth": 1, "alpha": 1, "color": "springgreen",
        #                        'edgecolor':'navy'})

        sns.histplot(area_values / 100, bins=num_bin, stat='count',
                     legend=False,
                     edgecolor='darkgreen',
                     facecolor='springgreen',
                     ax=axis_func)


        # Set y label
        axis_func.set_ylabel("Number of simulations", fontsize=15, labelpad=38)
        

    else:
        # Get another y axis - for 'Density'
        axis_density = axis_func.twinx()

        # Frequency plot
        # sns.distplot(area_values/unit, bins=num_bin, kde=False, ax=axis_func,
        #              hist_kws={"histtype": "bar", "linewidth": 1, "alpha": 1, "color": "springgreen",
        #                        'edgecolor':'darkgreen'})

        sns.histplot(area_values / 100, bins=num_bin, stat='count',
                     legend=False,
                     edgecolor='darkgreen',
                     facecolor='springgreen',
                     ax=axis_func)

        # Density plot
        sns.kdeplot(area_values / 100, ax=axis_density, color="r", linewidth=2)

        # Set y label for 'Frequency'
        axis_func.set_ylabel("Number of simulations", fontsize=20, labelpad=38)

        # Set y label for 'Density'
        axis_density.set_ylabel("Probability density", rotation=270, fontsize=20, labelpad=38)

        # Design ticks and x figures for 'Density'
        axis_density.tick_params(direction='out', length=8, pad=10)
        for item in (axis_density.get_xticklabels() + axis_density.get_yticklabels()):  # For x, y ticks' labels
            item.set_fontsize(15)

        # Approximately align with the first y axis. Please visit here for more information
        # https://stackoverflow.com/questions/45037386/trouble-aligning-ticks-for-matplotlib-twinx-axes
        # https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib (for
        # changing interval in ylim
        # Get lims of first axis - Frequency
        # len_axis1 = axis_func.get_ylim()

        start, end = axis_func.get_ylim()
        axis_func.yaxis.set_ticks(np.arange(start, end, 1))

        len_axis1 = axis_func.get_ylim()

        # Get lims of second axis - Density
        len_axis2 = axis_density.get_ylim()

        # Develop a function to get the general ticks for both x and y
        f = lambda x: len_axis2[0] + (x - len_axis1[0]) / (len_axis1[1] - len_axis1[0]) * (len_axis2[1] - len_axis2[0])

        # Get number of ticks
        num_ticks = f(axis_func.get_yticks())

        # Set the ticks for second y axis
        axis_density.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(num_ticks))

    # Area title
    title = f"Distribution of area values,\nresolution = {resolution_func} meters"
    title += add_title

    # Set up for x axis to make sure it appears with integers
    axis_func.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    # Ticks, title, and y label
    axis_func.tick_params(direction='out', length=8, pad=10)
    x_label = f'Areas (x100 '
    x_label += r'$\mathrm{m}^2$)'
    axis_func.set_xlabel(x_label, fontsize=20, labelpad=38)
    axis_func.set_title(title, fontsize=25, pad=25,
                        fontweight='bold')

    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)

    # Create text string
    textstr = "mean = {0:.3f}".format(area_dataframe_func.mean(axis=1)[0])
    textstr += "\nstdev = {0:.3f}".format(area_dataframe_func.std(axis=1)[0])

    # place a text box in upper left in axes coords
    # Refer here for more information
    # https://stackoverflow.com/questions/50869424/increase-line-separation-in-matplotlib-annotation-text/50888491
    axis_func.text(text_box_loc[0], text_box_loc[1], textstr, transform=axis_func.transAxes, fontsize=13,
                   linespacing=1.5,
                   fontweight= 'normal',
                   fontstyle='italic',
                   horizontalalignment='left',
                   verticalalignment='top')

    # Get max and min
    print("Maximum area:", np.max(area_values))
    print("Minimum area:", np.min(area_values))
    print("Mean of number of buildings:", area_dataframe_func.mean(axis=1)[0])
    print("Standard deviation of number of buildings:", area_dataframe_func.std(axis=1)[0])


def plot_building(axis_func, building_dataframe_func, resolution_func, text_box_loc,
                  density='y', add_title=""):
    """This function is to plot histogram of areas of simulations

    -----------
    References:
                https://stackoverflow.com/questions/69524514/how-to-modify-the-kernel-density-estimate-line-in-a-sns-histplot
                https://seaborn.pydata.org/generated/seaborn.distplot.html

                https://stackoverflow.com/questions/65400669/how-to-generate-two-separate-y-axes-for-a-histogram-on-the-same-figure-in-seabor
                https://stackoverflow.com/questions/26752464/how-do-i-align-gridlines-for-two-y-axis-scales-using-matplotlib
                https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.hist.html#matplotlib.axes.Axes.hist
                https://seaborn.pydata.org/generated/seaborn.kdeplot.html

                https://seaborn.pydata.org/generated/seaborn.histplot.html#seaborn.histplot
                https://seaborn.pydata.org/generated/seaborn.distplot.html
                https://stackoverflow.com/questions/27671748/how-to-print-y-axis-label-horizontally-in-a-matplotlib-pylab-chart
                https://stackoverflow.com/questions/24391892/printing-subscript-in-python

                https://stackoverflow.com/questions/45037386/trouble-aligning-ticks-for-matplotlib-twinx-axes (best
                answer for align two axis)
                https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib
    -----------

    -----------
    Arguments:
                axis_func:
                (axis in matplotlib)
                                        The subplot
                building_dataframe_func:
                (pandas dataframe)
                                        Dataframe of simulations' buildings
                resolution_func:
                (int or float)
                                        Resolution value in meter
                text_box_loc:
                (list)
                                        A list of text box coordinates [x, y]
                add_title:
                (string)
                                        Add more words into title
                density:
                (string)
                                        Performing "density" of "frequency"
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Get buildings' values under array format
    building_values = building_dataframe_func.iloc[0].to_numpy()

    # Number of bins/simulations
    num_bin = len(building_values)

    # Get histogram and frequency
    if density == 'y':
        # Draw kde (kernel density estimate) plot
        sns.kdeplot(building_values, color='crimson', ax=axis_func)

        # Set y label
        axis_func.set_ylabel("Probability density", fontsize=15, labelpad=38)


    elif density == 'n':
        # Draw histogram
        # sns.distplot(building_values, bins=num_bin, kde=False, ax=axis_func,
        #              hist_kws={"histtype": "bar", "linewidth": 1, "alpha": 1, "color": "deepskyblue", 'edgecolor':
        #                  'navy'})

        sns.histplot(building_values, bins=num_bin, stat='count',
                     legend=False,
                     edgecolor='navy',
                     facecolor='deepskyblue',
                     ax=axis_func)

        # Set y label
        axis_func.set_ylabel("Number of simulations", fontsize=15, labelpad=38)


    else:
        # Get another y axis - for 'Density'
        axis_density = axis_func.twinx()

        # Frequency plot
        # sns.distplot(building_values, bins=num_bin, kde=False, ax=axis_func,
        #              hist_kws={"histtype": "bar", "linewidth": 1, "alpha": 1, "color": "deepskyblue", 'edgecolor':
        #                  'navy'})

        sns.histplot(building_values, bins=num_bin, stat='count',
                     legend=False,
                     edgecolor='navy',
                     facecolor='deepskyblue',
                     ax=axis_func)

        # Density plot
        sns.kdeplot(building_values, ax=axis_density, color="r", linewidth=2)

        # Set y label for 'Frequency'
        axis_func.set_ylabel("Number of simulations", fontsize=20, labelpad=38)

        # Set y label for 'Density'
        axis_density.set_ylabel("Probability density", rotation=270, fontsize=20, labelpad=38)

        # Design ticks and x figures for 'Density'
        axis_density.tick_params(direction='out', length=8, pad=10)
        for item in (axis_density.get_xticklabels() + axis_density.get_yticklabels()):  # For x, y ticks' labels
            item.set_fontsize(15)

        # Approximately align with the first y axis. Please visit here for more information
        # https://stackoverflow.com/questions/45037386/trouble-aligning-ticks-for-matplotlib-twinx-axes
        # https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib (for
        # changing interval in ylim
        # Get lims of first axis - Frequency
        # len_axis1 = axis_func.get_ylim()

        start, end = axis_func.get_ylim()
        axis_func.yaxis.set_ticks(np.arange(start, end, 5))

        len_axis1 = axis_func.get_ylim()

        # Get lims of second axis - Density
        len_axis2 = axis_density.get_ylim()

        # Develop a function to get the general ticks for both x and y
        f = lambda x: len_axis2[0] + (x - len_axis1[0]) / (len_axis1[1] - len_axis1[0]) * (len_axis2[1] - len_axis2[0])

        # Get number of ticks
        num_ticks = f(axis_func.get_yticks())

        # Set the ticks for second y axis
        axis_density.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(num_ticks))

    # Set up for x axis to make sure it appears with integers
    axis_func.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    # Area title
    title = f"Distribution of number of buildings being inundated,\nresolution = {resolution_func} meters"
    title += add_title

    # Ticks, title, and y label
    axis_func.tick_params(direction='out', length=8, pad=10)
    axis_func.set_xlabel(fr'Number of buildings', fontsize=20, labelpad=38)
    axis_func.set_title(title, fontsize=25, pad=25,
                        fontweight='bold')

    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)

    # Create text string
    textstr = "mean = {0:.3f}".format(building_dataframe_func.mean(axis=1)[0])
    textstr += "\nstdev = {0:.3f}".format(building_dataframe_func.std(axis=1)[0])

    # place a text box in upper left in axes coords
    axis_func.text(text_box_loc[0], text_box_loc[1], textstr, transform=axis_func.transAxes, fontsize=13,
                   linespacing=1.5,
                   fontweight='normal',
                   fontstyle='italic',
                   horizontalalignment='left',
                   verticalalignment='top')

    # Get max and min
    print("Maximum number of buildings:", np.max(building_values))
    print("Minimum number of buildings:", np.min(building_values))
    print("Mean of number of buildings:", building_dataframe_func.mean(axis=1)[0])
    print("Standard deviation of number of buildings:", building_dataframe_func.std(axis=1)[0])


def find_nth(strings, value, n):
    """This function is to find the nth occurrence of the value in the strings

    -----------
    References: https://stackoverflow.com/questions/1883980/find-the-nth-occurrence-of-substring-in-a-string
    -----------

    -----------
    Arguments:
                strings:
                (string)
                            Strings is the where the value needs finding
                value:
                (int)
                            The value needs finding
                n:
                (int)
                            The ordinal number of occurrence of n
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Get the first occurrence of the value
    start = strings.find(value)

    # Start the loop to find the nth occurrence of the value
    while start >= 0 and n > 1:
        start = strings.find(value, start + len(value))
        n -= 1

    return start


def scatter_area_transformation(transformation_selection,
                                resolution_func,
                                axis_func, area_dataframe_func, add_title=""):
    """This function is to draw scatter plot between transformation and areas

    -----------
    References: https://www.digitalocean.com/community/tutorials/how-to-index-and-slice-strings-in-python-3
                https://stackoverflow.com/questions/31186019/rotate-tick-labels-in-subplot-pyplot-matplotlib-gridspec/52461208
                https://www.delftstack.com/howto/matplotlib/python-matplotlib-plot-superscript/
                https://stackoverflow.com/questions/9764298/how-to-sort-two-lists-which-reference-each-other-in-the-exact-same-way
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "tx" means x translation
                                            "ty" means y translation
                                            "c"  means combination
                resolution_func:
                (int or float)
                                            Resolution value in meter
                axis_func:
                (axis in matplotlib)
                                            The ordinal number of plot in subplot
                area_dataframe_func:
                (pandas dataframe)
                                            Dataframe of simulations' areas
                add_title:
                (string)
                                            Add more words into title
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Get list of areas
    areas_func = area_dataframe_func.iloc[0].to_numpy().tolist()

    # Get list of rotation
    if transformation_selection == "r":
        rotation_list = [float(angle[(find_nth(angle, "_", 1) + 1): find_nth(angle, "_", 2)]) for angle in
                         area_dataframe_func.columns]
        return_list = rotation_list

        # Set information for scatter plot
        title = f"Areas when transforming the grid - rotation,\nresolution = {resolution_func} meters"
        title += add_title

        # x label
        x_label = "Transformed simulations - rotation (degree)"


    # Get list of translation x
    elif transformation_selection == "tx":
        x_translation_list = [
            float(x_translation[(find_nth(x_translation, "_", 3) + 1): find_nth(x_translation, "_", 4)])
            for x_translation in area_dataframe_func.columns
        ]
        return_list = x_translation_list

        # Set information for scatter plot
        title = f"Areas when transforming the grid - east translation,\nresolution = {resolution_func} meters"
        title += add_title

        # x label
        x_label = "Transformed simulations - east translation (m)"

    # Get list of translation y
    elif transformation_selection == "ty":
        y_translation_list = [
            float(y_translation[(find_nth(y_translation, "_", 5) + 1):])
            for y_translation in area_dataframe_func.columns
        ]
        return_list = y_translation_list

        # Set information for scatter plot
        title = f"Areas when transforming the grid - north translation,\nresolution = {resolution_func} meters"
        title += add_title

        # x label
        x_label = "Transformed simulations - north translation (m)"

    # Get list of combination
    else:
        combination_list = [combination.replace("_", ":").replace("angle", "a") for combination in
                            area_dataframe_func.columns]
        return_list = combination_list

        # Set information for scatter plot
        title = f"Areas when transforming the grid - combination,\nresolution = {resolution_func} meters"
        title += add_title

        # x label
        x_label = "Transformed simulations - combination"

    # Get unit of area
    unit = np.power(resolution_func, 2)

    # Y labels
    y_label = f'Areas (x{unit} '
    y_label += r'$\mathrm{m}^2$)'

    # Sorting list
    new_return_list, new_areas_func = (list(val) for val in zip(*sorted(zip(return_list, areas_func))))

    # Convert list into string for x label
    string_return_list = [str(num) for num in new_return_list]

    # Draw scatter plot
    axis_func.scatter(string_return_list, new_areas_func/unit, s=250, facecolor="aqua", edgecolors="navy")

    # Set title
    axis_func.set_title(title, pad=25, fontsize=25, fontweight='bold')
    # Set x label
    axis_func.set_xlabel(x_label, fontsize=20, labelpad=38)
    # Set y label and y ticks
    axis_func.set_ylabel(y_label, rotation=-270, fontsize=20, labelpad=38)
    axis_func.tick_params('y', length=8, pad=10)
    # x label presentation
    axis_func.tick_params('x', labelrotation=90, length=8, pad=10)

    # Add grid
    axis_func.grid(which='both', axis='x', linestyle='--')

    # Control the grid
    for index, label in enumerate(axis_func.xaxis.get_ticklabels()):
        if index % 1 != 0:
            label.set_visible(False)

    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(15)

# END AREA PLOT ########################################################################################################



