# Prepare packages -----------------------------------------------------------------------------------------------------
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


def plotting_map(filtered_data_func, calculation_option,
                 axis_func, fig_func,
                 hex_list_func, colorbar_position,
                 extend_colorbar=None):
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
                extend_colorbar:
                                            Whether having or not having arrow to indicate
                                            the larger or smaller values in colorbar
                colorbar_position:
                                            "horizontal" or "vertical"
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Titles and labels:
    if calculation_option != 'cell':
        # Title for colorbar
        name_map = f"{calculation_option.capitalize()} of water depth (m)"

        # Title for contour map
        axis_func.set_title(f"{name_map} of water depth", pad=25, fontsize=25, fontweight='bold')
        axis_func.set_xlabel("NZTM, east (m)", fontsize=20, labelpad=38)
        axis_func.set_ylabel("NZTM, north (m)", fontsize=20, labelpad=38, rotation=-270)

        # Level information
        min_map_level = round_down(
            np.min(filtered_data_func[f'{calculation_option}'][filtered_data_func[f'{calculation_option}'] != -999]), 1)
        max_map_level = np.max(
            filtered_data_func[f'{calculation_option}'][filtered_data_func[f'{calculation_option}'] != -999] + 1)
        map_level = np.arange(min_map_level, max_map_level, 0.1)

    else:

        # Title for colorbar
        name_map = "Probability of each location being inundated (%)"

        # Title for contour map
        axis_func.set_title("Probability of each location being inundated", pad=25, fontsize=25, fontweight='bold')
        axis_func.set_xlabel("NZTM, east (m)", fontsize=20, labelpad=38)
        axis_func.set_ylabel("NZTM, north (m)", fontsize=20, labelpad=38, rotation=-270)

        # Max level
        min_map_level = 0
        max_map_level = 102
        map_level = np.arange(min_map_level, max_map_level, 1)

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
        map_colorbar.set_label(name_map, labelpad=30, fontsize=20)

    else:
        # Set up colorbar position
        cax = fig_func.add_axes([axis_func.get_position().x1 + 0.02,
                                axis_func.get_position().y0, 0.02,
                                axis_func.get_position().height])

        map_colorbar = fig_func.colorbar(contour_map_func, cax=cax, pad=5, ax=axis_func)
        map_colorbar.ax.tick_params(axis='y', direction='out', length=6, pad=10, labelsize=15)

        # Set up colorbar title
        map_colorbar.set_label(name_map, rotation=270, labelpad=38, fontsize=20)



def plotting_histogram(filtered_data_func, calculation_option,
                       axis_func, hex_list_func):
    """This function is to plot histogram of information regarding water depth

    -----------
    References: https://stackoverflow.com/questions/23061657/plot-histogram-with-colors-taken-from-colormap
                https://stdworkflow.com/67/attributeerror-rectangle-object-has-no-property-normed-solution
    -----------

    -----------
    Arguments:
                filtered_data_func:
                (pandas dataframe)
                                        Filtered dataset contains nodata values as -999
                calculation_option:
                (string)
                                        "mean" means to calculate mean
                                        "sd" means to calculate standard deviation
                                        "cv" means to calculate coefficient of variation
                                        "cell" means to calculate probability of each pixel being inundated
                axis_func:
                (axis in matplotlib)
                                        Subplot ordinal number of the contour map on the big plot
                hex_list_func:
                (string)
                                        A list of hex code colors (including or not including #)
    -----------

    -----------
    Returns:
                None.
    -----------

    """

    # Get z values
    z_func = filtered_data_func[f'{calculation_option}']

    # Select colormap
    hist_cm = get_gradient_cmap(hex_list_func)

    # Get the histogram
    hist_y_axis, hist_x_axis = np.histogram(z_func[z_func != -999], 50, density=False)

    hist_x_span = hist_x_axis.max() - hist_x_axis.min()

    hist_color = [hist_cm(((hist_x_val - hist_x_axis.min()) / hist_x_span)) for hist_x_val in hist_x_axis]

    axis_func.bar(hist_x_axis[:-1],
             hist_y_axis,
             color=hist_color,
             width=hist_x_axis[1] - hist_x_axis[0])

    for hist_item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):
        hist_item.set_fontsize(15)

    # Remove grid background lines (including x, y lines)
    axis_func.grid(False)
    axis_func.spines['top'].set_visible(False)
    axis_func.spines['right'].set_visible(False)
    axis_func.spines['bottom'].set_visible(False)
    axis_func.spines['left'].set_visible(False)

    # Set titles and labels
    if calculation_option != 'cell':
        # Title for colorbar
        name_map = f"{calculation_option.capitalize()} of water depth (m)"

        # Title for contour map
        axis_func.tick_params(direction='out', length=8, pad=10)
        axis_func.set_title(f"Histogram of {name_map} of water depth", pad=25, fontsize=25, fontweight='bold')
        axis_func.set_xlabel('Mean of depth values (m)', fontsize=20, labelpad=38)
        axis_func.set_ylabel('Frequency', fontsize=20, labelpad=38)

    else:

        # Title for contour map
        axis_func.tick_params(direction='out', length=8, pad=10)
        axis_func.set_title("Histogram of probability of each location being inundated", pad=25, fontsize=25,
                       fontweight='bold')
        axis_func.set_xlabel("Probability of each location being inundated (%)", fontsize=20, labelpad=38)
        axis_func.set_ylabel("Frequency", fontsize=20, labelpad=38)

# END STATISTICAL PLOT #################################################################################################



# AREA PLOT ############################################################################################################
def plot_area(axis_func, area_dataframe_func, density='y'):
    """This function is to plot histogram of areas of simulations

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                axis_func:
                (axis in matplotlib)
                                        The subplot
                area_dataframe_func:
                (pandas dataframe)
                                        Dataframe of simulations' areas
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

    # Number of bins/simulations
    num_bin = len(area_values)

    # Get histogram and frequency
    if density == 'y':
        sns.distplot(area_values, bins=num_bin, kde=True, ax=axis_func, hist=False)
    elif density == 'n':
        sns.distplot(area_values, bins=num_bin, kde=False, ax=axis_func)
    else:
        sns.distplot(area_values, bins=num_bin, kde=True, ax=axis_func)

    axis_func.set_xlabel(r'Areas ($m^{2}$)', fontsize=15, labelpad=30)
    axis_func.set_ylabel("Density", fontsize=15, labelpad=30)
    axis_func.set_title("Distribution of area values", fontsize=15, pad=30, fontweight='bold')

    # Get max and min
    print("Maximum area:", np.max(area_values))
    print("Minimum area:", np.min(area_values))


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


def scatter_area_transformation(transformation_selection, axis_func, area_dataframe_func):
    """This function is to draw scatter plot between transformation and areas

    -----------
    References: https://www.digitalocean.com/community/tutorials/how-to-index-and-slice-strings-in-python-3
                https://stackoverflow.com/questions/31186019/rotate-tick-labels-in-subplot-pyplot-matplotlib-gridspec/52461208
                https://www.delftstack.com/howto/matplotlib/python-matplotlib-plot-superscript/
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "tx" means x translation
                                            "ty" means y translation
                                            "c"  means combination
                axis_func:
                (axis in matplotlib)
                                            The ordinal number of plot in subplot
                area_dataframe_func:
                (pandas dataframe)
                                            Dataframe of simulations' areas
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
        rotation_list = [int(angle[(find_nth(angle, "_", 1) + 1): find_nth(angle, "_", 2)]) for angle in
                         area_dataframe_func.columns]
        return_list = rotation_list

        # Set information for scatter plot
        title = "Areas when transforming the grid - rotation"
        x_label = "Rotated angles (degree)"
        y_label = r'Areas ($m^{2}$)'


    # Get list of translation x
    elif transformation_selection == "tx":
        x_translation_list = [
            int(x_translation[(find_nth(x_translation, "_", 3) + 1): find_nth(x_translation, "_", 4)])
            for x_translation in area_dataframe_func.columns
        ]
        return_list = x_translation_list

        # Set information for scatter plot
        title = "Areas when transforming the grid - East translation"
        x_label = "East translations (m)"
        y_label = r'Areas ($m^{2}$)'


    # Get list of translation y
    elif transformation_selection == "ty":
        y_translation_list = [
            int(y_translation[(find_nth(y_translation, "_", 5) + 1):])
            for y_translation in area_dataframe_func.columns
        ]
        return_list = y_translation_list

        # Set information for scatter plot
        title = "Areas when transforming the grid - North translation"
        x_label = "East translations (m)"
        y_label = r'Areas ($m^{2}$)'

    # Get list of combination
    else:
        combination_list = [combination.replace("_", ":").replace("angle", "a") for combination in
                            area_dataframe_func.columns]
        return_list = combination_list

        # Set information for scatter plot
        title = "Areas when transforming the grid"
        x_label = "Combined transformation (m)"
        y_label = r'Areas ($m^{2}$)'

    # Draw scatter plot
    axis_func.scatter(return_list, areas_func)

    # Set title
    axis_func.set_title(title, pad=25, fontsize=15, fontweight='bold')
    # Set x label
    axis_func.set_xlabel(x_label, fontsize=15, labelpad=30)
    # Set y label
    axis_func.set_ylabel(y_label, rotation=-270, fontsize=15, labelpad=25)

    # Add grid
    axis_func.grid(which='both', axis='x', linestyle='--')

    # x label presentation of combined transformation
    if transformation_selection == "c":
        axis_func.tick_params('x', labelrotation=90)

    # Control the grid
    for index, label in enumerate(axis_func.xaxis.get_ticklabels()):
        if index % 3 != 0:
            label.set_visible(False)

    for item in (axis_func.get_xticklabels() + axis_func.get_yticklabels()):  # For x, y ticks' labels
        item.set_fontsize(13)

# END AREA PLOT ########################################################################################################



