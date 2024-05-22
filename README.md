# Grid_Orientation
Quantifying uncertainty in flood predictions caused by transforming square grid orientation in LiDAR-derived DEM process
Authors: M. Nguyen, M. D. Wilson, E. M. Lane, J. Brasington, and R. A. Pearson

## Summary

Digital elevation datas, or DEMs, are critical for generating reliable flood predictions. The most common method to generate DEMs involves sampling and interpolating LiDAR data onto a North-South square grid. However, the orientation of this grid, which can introduce variability in elevation data and thus influence flood predictions, is frequently overlooked. The below figure shows the difference in sampling an artificial LiDAR point cloud with all ground elevation values of 1 m and river elevation values of 0.5 m onto a North-South square grid versus a 45-degree rotated grid using 10-meter resolution. This resulted in significant differences in the representation of the river, as seen in the river profiles, which can create different waterways along its length. 

![problem_idea_003](https://github.com/Martin20494/Grid_Orientation/assets/55137629/65c5d839-0db1-4d79-aecd-7baa78c5b4a5)

Due to its potential significant impacts on the flood modelling process, in this study, we quantify the uncertainty in flood model outputs caused by transforming this grid orientation. Using a Monte Carlo method, we produced multiple DEMs by randomly rotating and/or translating the square grid orientation to predict floods for uncertainty analysis. This Monte Carlo framework was also applied at various resolutions (2, 5, 10, and 20 meters) and flood return periods (5, 10, 20, 50, and 1000 years). 

Results show that the highest uncertainty in flood predictions occurred predominantly at the flood extent boundaries and near the river. Rotating the grid orientation caused more uncertainty than translating. The highest uncertainty was found when both rotating and translating the grid. The variability in both flooded areas and the number of flooded buildings followed these patterns. Finer resolutions revealed fewer variations in flood predictions and less variability in flooded areas and the number of flooded buildings. These variations were also determined by the amount of river discharge. Significant variations occurred depending on whether the river discharge was insufficient to cover the floodplain surfaces. Also, the flooded areas and the number of flooded buildings could expand or contract based on whether the water reached new locations and these places included residential buildings.


![S3_proportion_wd](https://github.com/Martin20494/Grid_Orientation/assets/55137629/840e5b4e-4801-43e2-80ff-0804d696ccbb)

The above figure shows the proportions of times a location was predicted to be flooded in all simulations. It usggests that in some simulations, topopographic effects, caused by the changes in the grid orientation, prevented water from reaching some specific locations, such as the upper part of the blue zoomed-in image in the figure.

## Workflow

![transformation_process](https://github.com/Martin20494/Grid_Orientation/assets/55137629/4ecdd3b5-2e28-41b2-ae65-8ba4044b20d8)

## Data

The DEMs in this work is created thanks to the GeoFabrics package and instruction provided by [R.A.Pearson](https://github.com/rosepearson/GeoFabrics/wiki) her work was published and can be read [here](https://www.sciencedirect.com/science/article/pii/S1364815223002281). Four files attached here were used to run the simulations.

#### 1. File validation_calibration

This file includes three other files - data, results, and scripts. These files were used to validate and calibrate the LISFLOOD-FP flood model.
- Data file includes all data that are necessary for validating and calibrating the flood model. Notice here the land (51559) and bathymetry (50554) were downloaded on 15/04/2024. The OpenStreetMap data was downloaded into cache file on 25/01/2023. Hence, the reproducible work could be different because the versions of these data could be changed.
- Scripts file includes all scripts to validate and calibrate the models.
- Results file includes the results from our work after running the scripts file.

#### 2. File simulation

This file includes two files - data and script.
- Data file includes initial data which cannot be generated automatically from the script file.
- Script file includes three other files - Modelling, Analysing, and Executing. The Modelling file includes all the modules used to create multiple DEMs and flood predictions by transforming the grid orientation. The Analysing file includes all the modules used to extract all uncertainties and analyse them. The Executing file includes all the scripts to run the Modelling and Analysing file.
- The results for these scripts are very large (4TB). Please contact with the author if you would like to know more details.

#### 3. File other_files

This file includes two big files - script and data.
- Data file includes all the files used for figures in the journal.
- Script file includes all the script used to create data file.

#### 4. File comparison

- Comparison_script file is used to run the codes using to compare all transformation types, resolutions, and flood return periods.
- Comparison_result files is the result after running the script.

## Acknowledgement
This study is part of the NIWA-led national flood hazard assessment programme, ["Reducing flood inundation hazard and risk across Aotearoa - New Zealand"](https://niwa.co.nz/hazards/ma-te-haumaru-o-nga-puna-wai-o-rakaihautu-ka-ora-mo-ake-tonu-increasing-flood), funded by the Ministry for Business, Innovation and Employment (MBIE) Endeavour Programme. A PhD scholarship is provided from this programme to M. Nguyen. The LISFLOOD-FP flood model was developed at the University of Bristol and is available [here](https://www.seamlesswave.com/LISFLOOD8.0)

## Contacts

Github author: Martin Nguyen martinnguyen20494@gmail.com



