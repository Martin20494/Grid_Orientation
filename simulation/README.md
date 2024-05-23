Download ```data``` and ```scripts``` folders into your local folder. 

#### 1. File data

The ```data``` file includes ```0_lidar_data``` and ```other_data```. These files should be put into a simulation version folder (e.g: Put them into folder version_001 used to generate simulations of rotation and North-East translation).

- ```0_lidar_data``` includes ```flow_and_friction.csv.gz``` and ```rec2_3.geojson``` files (river features) provided by NIWA to estimate the river bathymetry
- ```other_data``` includes shape file of observed flood debris (January-2005 flood event at Waikanae River) and linestring represents for tide boundary used for flood modelling.

#### 2. File scripts

The ```scripts``` file includes ```Analysing```, ```Modelling```, and ```Executing``` folders. These file are used for generating simulations. 

- ```Modelling``` and ```Analysing```: Includes two separate groups of modules to generate and then analyse the simulations. ```Modelling``` will create multiple DEMs from LiDAR data downloaded from the [OpenTopography] (https://portal.opentopography.org/datasets). These DEMs will then be used in LISFLOOD-FP model to produce multiple water depths and water surfaces. These data will then be analysed in ```Analysing```. Please change the ```versionModule``` to the name of the simulation version you are running. For example, I named ```version_001``` for rotation and North-Eastranslation simulation. Also change the ```MAIN_DIR``` to your own path.
- ```Executing``` folder: Includes ```Analysing```, ```Comparing```, and ```Modelling``` files used to run modules in ```Analysing``` and ```Modelling``` folders to generate and compare simulations. These files were written in jupyter notebook format. Therefore, we suggest to also use jupyter notebook with the environment installed earlier. Please follow the instructions for generating each simulation version for different transformation types, resolutions, and flood return periods written in the script. There are 13 version of these different scenrios applied on this uncertainty in total. Please change the ```os.chdir``` in each file into the folder containing ```Analysing``` or ```Modelling``` folders.

**_Notice_**: This project was run on CPU machine with 8 cores and 64GB memory, you might want to choose different ```size_of_processor``` and ```size_of_chunk``` in ```Modelling``` files according to your machine features. It will take 1-2 days to generate a version but some might take more than a week (e.g: 2-meter resolution). To make the work faster, the Github author used [Nesi computers](https://www.nesi.org.nz/services/high-performance-computing-and-data-analytics) to generate multiple simulations of flood predictions at the same time after creating DEMs, if you would like to proceed the same please visit [other_files](https://github.com/Martin20494/Grid_Orientation/tree/main/other_files) for more information. Please contact with the Github author if further information is needed. The image below shows how 13 versions and files in a version looks like after running these files.

<div align="center">
	<img width = "90%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/simulation/folder_example/folders_examples.jpg">
</div>
