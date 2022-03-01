# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                            # For paths of sub-folders
import shutil                                   # For file handling
import pathlib                                  # For creating and manipulating the directory path
# ----------------------------------------------------------------------------------------------------------------------


def save_plot(transformation_selection, fig_func, plot_name, extension, dpi_func, addition_path=""):
    """This function is to save the plots

    -----------
    References: https://stackoverflow.com/questions/4325733/save-a-subplot-in-matplotlib
                https://stackoverflow.com/questions/16032389/pad-inches-0-and-bbox-inches-tight-makes-the-plot-smaller-than-declared-figsiz
                https://stackoverflow.com/questions/10041627/how-to-make-savefig-save-image-for-maximized-window-instead-of-default-size
    -----------

    -----------
    Arguments:
                transformation_selection:
                (string)
                                            "r" means rotation
                                            "t" means translation
                                            "c" means combination
                fig_func:
                (figure in matplotlib)
                                            Figure in subplot matplotlib
                plot_name:
                (string)
                                            Name of the plot
                extension:
                (string)
                                            "png" or "pdf"
                dpi_func:
                (int or float)
                                            The resolution in dots per inch. Please visit here for more information
                                            'https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html'
                addition_path:
                (string)
                                            Path for saving plots
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    if transformation_selection == "r":
        saving_path = plot_rotation
    elif transformation_selection == "t":
        saving_path = plot_translation
    elif transformation_selection == "c":
        saving_path = plot_combination
    else:
        saving_path = addition_path

    fig_func.savefig(fr"{saving_path}\\{plot_name}.{extension}", bbox_inches='tight', dpi=dpi_func)


def graph_folder(oldpath,
                 plotname,
                 newpath,
                 newname):
    """This function is to copy files from a folder to another folder. Mainly used for graphs

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                oldpath:
                (string)
                                Old directory
                plotname:
                (string)
                                Name of plot including extension
                newpath:
                (string)
                                New directory
                newname:
                (string)
                                New name for the plot
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Create folder if not exist
    pathlib.Path(f"{newpath}").mkdir(parents=True, exist_ok=True)

    # Copy graph
    shutil.copy2(f"{oldpath}\\{plotname}.png", f"{newpath}\\{newname}.png")


def get_figure_folder(
        super_selection,
        selection="",
        version_func="",
        addition=""):

    """This function is to set up files for copying. Mainly used for graphs

    -----------
    References:
                None.
    -----------

    -----------
    Arguments:
                super_selection:
                (string)
                                    Options between 'no comparsion' and 'comparison'
                selection:
                (string)
                                    Options between 'transformation' and 'resolution'
                version_func:
                (string)
                                    Name of version that contain information which needs extracting
                addition:
                (string)
                                    Options between '' and ' histogram'
    -----------

    -----------
    Returns:
                None.
    -----------

    """
    # Set up drive and main_folder which contain all necessary information
    drive = "S"
    main_folder = "LISFLOOD"
    newname = ""

    # Choose 'no comparison' or 'comparison'
    if super_selection == 'no comparison':

        # Get information from 'mean', 'sd', 'cv', 'cell.
        # If 'histogram', 'area', 'building' will be added
        if addition == "":
            statistical_list_func = ['mean', 'sd', 'cv', 'cell']
        else:
            statistical_list_func = ['mean', 'sd', 'cv', 'cell', "area", "building"]

        # Loop to get information
        for newfolder in statistical_list_func:
            oldname = newfolder + addition

            # Choose between 'resolution' or 'transformation'
            if selection == 'resolution':
                if version_func == "version_21":
                    newname = f"2m_{oldname}"
                elif version_func == "version_26":
                    newname = f"5m_{oldname}"
                elif version_func == "version_12":
                    newname = f"10m_{oldname}"
                else:
                    newname = f"20m_{oldname}"

            elif selection == 'transformation':
                if version_func == "version_12":
                    newname = f"combination_{oldname}"
                elif version_func == "version_17":
                    newname = f"rotation_{oldname}"
                elif version_func == "version_7":
                    newname = f"translation_{oldname}"
                elif version_func == "version_8":
                    newname = f"xtranslation_{oldname}"
                else:
                    newname = f"ytranslation_{oldname}"

            header_1 = fr"{drive}:\\{main_folder}\\{version_func}"
            oldpath = f"{header_1}\\7_results\\un_combination\\uncom_plot"
            newpath = fr"P:\\Courses\\PhD\\PhD - Thesis information\\PhD documents\\Report\\LaTex\\figure\\{newfolder}"

            graph_folder(
                oldpath,
                oldname,
                newpath,
                newname)

    else:
        for version_selection in ['resolution', 'transformation']:
            for oldfolder in ['boxplot', 'density']:
                if oldfolder == "density":
                    statistical_list_func = ['mean', 'sd', 'cv', 'cell']
                else:
                    statistical_list_func = ['mean', 'sd', 'cv', 'cell', "area", "building"]
                for newfolder in statistical_list_func:
                    oldname = f"{newfolder}_{version_selection}"
                    newname = f"{newfolder}_{version_selection}_{oldfolder}"

                    header_2 = fr"{drive}:\\{main_folder}\\{oldfolder}"
                    oldpath = f"{header_2}\\{version_selection}"
                    newpath = fr"P:\\Courses\\PhD\\PhD - Thesis information\\PhD documents\\Report\\LaTex\\figure\\{newfolder}"

                    graph_folder(
                        oldpath,
                        oldname,
                        newpath,
                        newname
                    )