Download `data` and `scripts` folders into your local computer. 

#### 1. File data

The `data` folder: Includes `0_lidar_data` and `other_data` sub-folders. These files should be put into a simulation version folder (e.g: Put them into folder version_001 used to generate simulations of rotation and translation).

- `0_lidar_data`: Includes `flow_and_friction.csv.gz` and `rec2_3.geojson` files (river features) provided by NIWA to estimate the river bathymetry. These data can be found [here](https://data-niwa.opendata.arcgis.com/apps/NIWA::new-zealand-river-flood-statistics-app/explore).

- `other_data`: includes shape file of observed flood debris (January-2005 flood event at Waikanae River) and linestring represents for tide boundary used for flood modelling.

#### 2. File scripts

The `scripts` folder includes `Analysing`, `Modelling`, and `Executing` sub-folders. These sub-folders are used for generating simulations. 

- `Modelling` and `Analysing`: Includes two separate groups of modules to generate and then analyse the simulations. `Modelling` will create multiple DEMs from LiDAR data downloaded from the [OpenTopography](https://portal.opentopography.org/datasets). These DEMs will then be used in LISFLOOD-FP model to produce multiple water depths and water surfaces. These data will then be analysed in `Analysing`.

- `Executing`: Includes `Analysing`, `Comparing`, and `Modelling` files used to run modules in `Analysing` and `Modelling` folders to generate and compare simulations. These files were written in jupyter notebook format. Therefore, we suggest to also use jupyter notebook with the environment installed earlier. Please follow the instructions for generating each simulation version for different transformation types, resolutions, and flood return periods written in the script. Please change the `os.chdir` in each file into the folder containing `Analysing` or `Modelling` folders.

#### 3. Usage

After [installation](https://github.com/Martin20494/Grid_Orientation?tab=readme-ov-file#environment-installation), download these subfolders to a local folder. In `Modelling` and `Analysing`, change the `MAIN_DIR` (path to your local folder) in module `folder.py`, change the name of simulation version in module `versionModule`. For example, I named `version_001` for rotation and translation with 10-m grid and using January-2015 event. You might also want to change the [API key](https://www.linz.govt.nz/guidance/data-service/linz-data-service-guide/web-services/creating-api-key) to crawl data from LINZ in `dataPreparation.py` in `Modelling`.

There are 13 simulation versions, each of them represents for a set of 50 transformation values, a resolution, and a flood event. To create 50 simulations for each version, open `Executing` folder and access to jupyter notebook `Modelling.ipynb`, change things as following bullet points befor running all cells:

	- Change the path in `0. Change directory` for cell directing to where the `Modelling.ipynb` is stored in your local computer.
 
 	- Change the `2. Necessary variables` to the simulation version you want to run.

	- Change the `2.2. Other basic variables` (probably the second cell) and `2.4. Preparing some flood model inputs` to the size and chunk to suit your machine.

	- Change the `2.5. Random transformations` to the simulation version you want to run.

 	- Change the number of processors in all cells of `3. Execution` to suit your machine.

To generate analysis result, in `Executing` folder, access to jupyter notebook `Analysis.ipynb`, change things as following bullet points before running all cells:

	- Change the parth in `0. Change directory` for cell directing to where the `Analysis.ipynb` is stored in your local computer.

 	- Change the `2. Data preparation` (probably the resolution) to the simulation version you want to run.

  	- Change the number of processors in all cells of `3. Generate csv files` to suit your machine.

After generating all (13) subsets, to compare, in `Executing` folder, access to jupyter notebook `Comparison.ipynb`, change things as following bullet points before running all cells:

	- Change the parth in the very first cell to let it direct to where the `Analysis.ipynb` is stored in your local computer.

 	- Change the paths in variables `trans_list_filename`, `res_list_filename`, and `events_list_filename` into your local paths where you store all the subsets. Their orders should be according to the `trans_name`, `res_name`, and `event_name`.

  	- Change the paths to store all the results, they are in cells with the last line of codes with comment line `# Save fig`.

**_Notice_**: This project was run on CPU machine with 8 cores and 64GB memory, you might want to choose different `size_of_processor` and `size_of_chunk` in `Modelling` files according to your machine features. It will take 1-2 days to generate a version but some might take more than a week (e.g. 2-meter resolution). Please contact with the Github author if further information is needed. The image below shows how 13 versions and files in a version look like after running these files.

<div align="center">
	<img width = "90%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/simulation/folder_example/folders_examples.jpg">
</div>
