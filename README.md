# Grid Alignment

Title: Quantifying uncertainty in flood predictions LISFLOOD-FP due to arbitrary conventions in grid alignment in LiDAR-derived digital elevation models

Authors: Martin Nguyen, Matthew D. Wilson, Emily M. Lane, James Brasington, and Rose A. Pearson

## Introduction

Digital elevation models (DEMs) generated by sampling and interpolating LiDAR data onto a square grid can produce reliable flood predictions. However, the arbitrary conventions in grid alignment that can introduce uncertainty in flood predictions are frequently overlooked. To demonstrate this, we illustrates a simple case of an artificial LiDAR point cloud with all ground elevation values of 1 m and a 10 m wide river with elevation values of 0.9 m meandering through it. These LiDAR data were sampled and interpolated onto a square grid versus a versus a $45^{\circ}$ rotated square grid using 10-m spatial resolution. The differences resulted from these grids are demonstrated in the river profiles, as shown in the below figure. These differences may affect the depth and subsequently impact the floodplain once it floods. Hence, this project aims to quantify the uncertainty in flood predictions arising from arbitrary conventions in this square grid alignment.

<div align="center">
	<img width = "50%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data/Problem/problem_idea.jpg">
</div>

## Site study

The Waikanae river, situated on the West Coast of the Wellington Region in New Zealand, was selected for this project. Figure (a) shows the topographic features of the Waikanae River extending from the Waikanae Treatment Plant gauge approximately 7 km towards the coastal region. The DEM in Figure (b) and the roughness length in Figure (c) were generated by using GeoFabrics package [(version 0.9.4)](https://github.com/rosepearson/GeoFabrics/commit/7a4ebca4e3f2e0c3f9218f7b1ef07ee470d51cb7) provided by [Dr. Rose A. Pearson](https://github.com/rosepearson/GeoFabrics/wiki) and published [here](https://www.sciencedirect.com/science/article/pii/S1364815223002281)

<div align="center">
	<img width = "60%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data/DEM_hillshade/waikanae_all.jpg">
</div>

## Workflow

We used LISFLOOD-FP flood model within a Monte Carlo framework to generate multiple flood predictions based on the January-2005 Waikanae River flood event for analysis. This approach was applied to five transformation types (East translation, North translation, translation, rotation, and a combination of rotation and translation), four resolutions (2, 5, 10, and 20 meters), and six flood scenarios (with return periods of 5, 10, 20, 50, 80, and 1000 years).

<div align="center">
	<img width = "40%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data/Simulation_process/transformation_list.jpg">
	<img width = "70%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data/Simulation_process/transformation_process.jpg">
</div>

## Representative result

The figure below demonstrates the effect of transforming the square grid orientation on the flood predictions. It shows the proportion of times a location/each pixel was predicted to be flooded across all simulations. The results suggest that in some simulations, topographic effects due to the changes in the grid orientation prevented water from reaching specific locations. This is evident in the upper part of the blue zoomed-in image or at the flood extent boundaries in the figure.

<div align="center">
	<img width = "60%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data/All_results/Results_representative/S3_proportion_wd.jpg">
</div>

## Reproducibility

The project was written under modules format and has not been updated to class format yet. Therefore, the users might need to download all the modules to reproduce this project. After the instruction here, please click on each folder for further guidelines. If further information is needed please contact with the Github author.

### Environment installation

Download file ```grid_orientaion_packages.yaml``` and use the following command (reference from the answer of [@merv](https://stackoverflow.com/questions/76800978/conda-invalidversionspec-invalid-version-error-when-tryin-to-install-from-requi)) to recreate the anaconda environment (only support Windows at the moment).

```
conda env create -n gridorientation -f gridorientation.yaml
```

### Data

This project includes 3 main folders as described below.

- [calibration](https://github.com/Martin20494/Grid_Orientation/tree/main/validation_calibration): Including `data`, `scripts`, and `results` sub-folders used to calibrate the LISFLOOD-FP flood model. 

- [simulation](https://github.com/Martin20494/Grid_Orientation/tree/main/simulation): Including `data` and `scripts` sub-folders used to generate simulations.

- [other_files](https://github.com/Martin20494/Grid_Orientation/tree/main/other_files): Including `data` and `scripts` sub-folders used to store all results from this project, figures for manuscript, and data sources.

## Acknowledgement

This project is part of the NIWA-led national programme, ["Reducing flood inundation hazard and risk across Aotearoa - New Zealand"](https://niwa.co.nz/hazards/ma-te-haumaru-o-nga-puna-wai-o-rakaihautu-ka-ora-mo-ake-tonu-increasing-flood), funded by the Ministry for Business, Innovation and Employment (MBIE) Endeavour Programme. A PhD scholarship was provided from this programme to Martin Nguyen. The LISFLOOD-FP flood model was developed at the University of Bristol and is available [here](https://www.seamlesswave.com/LISFLOOD8.0).

## Contacts

Github author: Martin Nguyen martinnguyen20494@gmail.com



