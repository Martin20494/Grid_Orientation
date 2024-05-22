# Grid_Orientation

Title: Quantifying uncertainty in flood predictions caused by transforming square grid orientation in LiDAR-derived DEM process

Authors: Martin Nguyen, Matthew D. Wilson, Emily M. Lane, James Brasington, and Rose A. Pearson

## Introduction

Digital elevation models (DEM), typically generated by sampling and interpolation LiDAR data onto a North-South square grid, are crucial for reliable flood predictions. However, the orientation of this grid, which can introduce variability in elevation data and impact flood predictions, is often overlooked. The below figure illustrates the difference in sampling and interpolating artificial LiDAR data (ground elevation: 1 m, river elevation: 0.5 m) onto a North-South square grid versus a 45-degree rotated grid at 10-meter resolution. This significantly affects river representation as seen in river profiles, potentially creating varying waterways. Hence, this project aims to quantify the uncertainty in flood model outputs caused by transforming this grid orientation.

<div align="center">
	<img width = "65%" src="https://github.com/Martin20494/Grid_Orientation/assets/55137629/65c5d839-0db1-4d79-aecd-7baa78c5b4a5)](https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data_forpublication/Problem/problem_idea_003.jpg">
</div>

## Workflow

We used a traditional Monte Carlo framework to generate multiple DEMs by rotating and/or translating a square grid. These DEMs were created using GeoFabrics package [(version 0.9.4)](https://github.com/rosepearson/GeoFabrics/commit/7a4ebca4e3f2e0c3f9218f7b1ef07ee470d51cb7) provided by [Dr. Rose A. Pearson](https://github.com/rosepearson/GeoFabrics/wiki) and published [here](https://www.sciencedirect.com/science/article/pii/S1364815223002281). Using the LISFLOOD-FP flood model, we produced multiple flood predictions for analysis from these DEMs. This approach was applied to various transformations (East, North, North-East translation, rotation, and combined rotation and North-East translation), resolutions (2, 5, 10, and 20 meters), and flood return periods (5, 10, 20, 50, 80 (or January-2005), and 1000 years).

![transformation_process](https://github.com/Martin20494/Grid_Orientation/assets/55137629/4ecdd3b5-2e28-41b2-ae65-8ba4044b20d8)

## Representative result

The figure below demonstrates the effect of transforming the square grid orientation on the flood predictions. It shows the proportion of times a location/each pixel was predicted to be flooded across all simulations. The results suggest that in some simulations, topographic effects due to the changes in the grid orientation prevented water from reaching specific locations. This is evident in the upper part of the blue zoomed-in image or at the flood extent boundaries in the figure.

<div align="center">
	<img width = "90%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data/All_results/Results_representative/S3_proportion_wd.jpg">
</div>

## Reproducibility

The project was written under modules format and has not been updated to class format yet. Therefore, the users might need to download all the modules to reproduce this project. After the instruction here, please click on each folder for further guidelines. If further information is needed please contact with the Github author.

### Environment installation

Download file ```gri_orientaion_packages.yaml``` and use the following command (reference from the answer of [@merv](https://stackoverflow.com/questions/76800978/conda-invalidversionspec-invalid-version-error-when-tryin-to-install-from-requi)) to recreate the anaconda environment (only support Windows at the moment).

```
conda env create -n gridorientation -f gridorientation.yaml
```

### Data

This project includes 3 main folders as described below.

- [validation_calibration](https://github.com/Martin20494/Grid_Orientation/tree/main/validation_calibration): Including data, scripts, and results sub-folders used to validate and calibrate the LISFLOOD-FP flood model. 

- [simulation](https://github.com/Martin20494/Grid_Orientation/tree/main/simulation): Including data and scripts sub-folders used to generate simulations.

- [other_files](https://github.com/Martin20494/Grid_Orientation/tree/main/other_files): Including data and scripts sub-folders used to store all results from this project, figures for manuscript, and data sources.

## Acknowledgement
This project is part of the NIWA-led national programme, ["Reducing flood inundation hazard and risk across Aotearoa - New Zealand"](https://niwa.co.nz/hazards/ma-te-haumaru-o-nga-puna-wai-o-rakaihautu-ka-ora-mo-ake-tonu-increasing-flood), funded by the Ministry for Business, Innovation and Employment (MBIE) Endeavour Programme. A PhD scholarship was provided from this programme to Martin Nguyen. The LISFLOOD-FP flood model was developed at the University of Bristol and is available [here](https://www.seamlesswave.com/LISFLOOD8.0).

## Contacts

Github author: Martin Nguyen martinnguyen20494@gmail.com



