# Prepare packages -----------------------------------------------------------------------------------------------------
import numpy as np                                                   # For data manipulation
import matplotlib.colors as mcolors                                  # For coloring
# ----------------------------------------------------------------------------------------------------------------------


def hex_to_rgb(value_func):
    """This function is to convert hex to red green blue (rgb) colors.
    The code is based on the link in reference (by Kerry Halupka)

    -----------
    References: https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
    -----------

    -----------
    Arguments:
                value_func:
                (string)
                                 Hex color code under 6-characters string format
    -----------

    -----------
    Returns:
                rgb_value_func:
                (tuple)
                                 A tuple of rgb values with length of 3
    -----------

    """
    # Remove '#' in string
    hex_value_func = value_func.strip('#')

    # Get the quantity of colors
    level_func = len(hex_value_func)

    # Convert color from hex to rgb and store in a tuple
    rgb_value_func = tuple(
        int(hex_value_func[i:i + level_func // 3], 16) for i in range(0, level_func, level_func // 3))

    return rgb_value_func



def rgb_to_dec(value_func):
    """This function is to convert rgb to decimal colors (dividing each value by 256).
    The code is based on the link in reference (by Kerry Halupka)

    -----------
    References: https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
    -----------

    -----------
    Arguments:
                value_func:
                (tuple)
                                A tuple of rgb color code (from 0 to 256) with length of 3
    -----------

    -----------
    Returns:
                dec_value_func:
                (list)
                                A list of color decimal values with length of 3
    -----------

    """
    return [val / 256 for val in value_func]



def get_gradient_cmap(hex_list_func, float_list_func=None):
    """This function is to create gradient colors. The code is based on the link in reference (by Kerry Halupka)
    If float_list_func is None, colour map will graduate linearly between each color in hex_list
    If float_list_func is not None, each color in hex_list_func is mapped to the respective location in float_list_func

    -----------
    References: https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
    -----------

    -----------
    Arguments:
                hex_list_func:
                (list)
                                                A list of hex code color under string format
                float_list_func:
                (list)
                                                A list of floats (between 0 and 1), same length as hex_list_func.
                                                Must start with 0 and end with 1.
    -----------

    -----------
    Returns:
                colour_map:
                (color map in matplotlib)
                                                Color under matplotlib color map
    -----------

    """
    # Get rgb list
    rgb_list = [rgb_to_dec(hex_to_rgb(color_code)) for color_code in hex_list_func]

    # Check float list
    if float_list_func:
        pass
    else:
        float_list_func = list(np.linspace(0, 1, len(rgb_list)))

    # Build up gradient colors
    color_dict = dict()
    for number, color in enumerate(['red', 'green', 'blue']):
        color_list = [[float_list_func[i], rgb_list[i][number], rgb_list[i][number]] for i in
                      range(len(float_list_func))]
        color_dict[color] = color_list

    color_map = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=color_dict, N=256)

    return color_map
