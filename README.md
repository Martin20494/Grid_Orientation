# Grid_Orientation

Title: Quantifying uncertainty in flood predictions caused by transforming square grid orientation in LiDAR-derived DEM process

Authors: Martin Nguyen, Matthew D. Wilson, Emily M. Lane, James Brasington, and Rose A. Pearson

## Introduction

Digital elevation models (DEM), typically generated by sampling and interpolation LiDAR data onto a North-South square grid, are crucial for reliable flood predictions. However, the orientation of this grid, which can introduce variability in elevation data and impact flood predictions, is often overlooked. The below figure illustrates the difference in sampling and interpolating artificial LiDAR data (ground elevation: 1 m, river elevation: 0.5 m) onto a North-South square grid versus a 45-degree rotated grid at 10-meter resolution. This significantly affects river representation as seen in river profiles, potentially creating varying waterways. Hence, this project aims to quantify the uncertainty in flood model outputs caused by transforming this grid orientation.

<div align="center">
	<img width = "65%" src="https://github.com/Martin20494/Grid_Orientation/assets/55137629/65c5d839-0db1-4d79-aecd-7baa78c5b4a5)](https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data_forpublication/Problem/problem_idea_003.jpg">
</div>

## Workflow

We used the traditional Monte Carlo framework to generate multiple DEMs by rotating and/or translating a square grid. Using the LISFLOOD-FP flood model, we produced multiple flood predictions for analysis. This approach was applied to various transformations (East, North, North-East translation, rotation, and combined rotation and North-East translation), resolutions (2, 5, 10, and 20 meters), and flood return periods (5, 10, 20, 50, 80 (or January-2005), and 1000 years).

![transformation_process](https://github.com/Martin20494/Grid_Orientation/assets/55137629/4ecdd3b5-2e28-41b2-ae65-8ba4044b20d8)

## Representative result

The figure below demonstrates the effect of transforming the square grid orientation on the flood predictions. It shows the proportion of times a location/each pixel was predicted to be flooded across all simulations. The results suggest that in some simulations, topographic effects due to the changes in the grid orientation prevented water from reaching specific locations. This is evident in the upper part of the blue zoomed-in image or at the flood extent boundaries in the figure.

<div align="center">
	<img width = "90%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data/All_results/Results_representative/S3_proportion_wd.jpg">
</div>

## Reproducible work

The instruction below is for reproducing this project. It will take 1-2 days to generate a version but some might take more than a week (e.g: 2-meter resolution). Please contact with the Github author if further information is needed.

### Environment installation

Download file ```gri_orientaion_packages.yaml``` and use the following command (reference from the answer of [@merv](https://stackoverflow.com/questions/76800978/conda-invalidversionspec-invalid-version-error-when-tryin-to-install-from-requi)) to recreate the anaconda environment (only support Windows at the moment).

```
conda env create -n gridorientation -f gridorientation.yaml
```

### Data

The DEMs in this work are created thanks to the GeoFabrics package (version 0.9.4) provided by [Dr. Rose A. Pearson](https://github.com/rosepearson/GeoFabrics/wiki) her work can be read [here](https://www.sciencedirect.com/science/article/pii/S1364815223002281). This project includes 3 main folders as described below.

#### 1. File validation_calibration

Data, scripts, and results sub-folders in this folder were used to validate and calibrate the LISFLOOD-FP flood model. 

#### 2. File simulation

Data and scripts sub-folders in this folder were used to generate simulations.

#### 3. File other_files

Data and scripts sub-folders in this folder were used to store all results from this project, figures for manuscript, and data sources.

## Acknowledgement
This project is part of the NIWA-led national flood hazard assessment programme, ["Reducing flood inundation hazard and risk across Aotearoa - New Zealand"](https://niwa.co.nz/hazards/ma-te-haumaru-o-nga-puna-wai-o-rakaihautu-ka-ora-mo-ake-tonu-increasing-flood), funded by the Ministry for Business, Innovation and Employment (MBIE) Endeavour Programme. A PhD scholarship was provided from this programme to Martin Nguyen. The LISFLOOD-FP flood model was developed at the University of Bristol and is available [here](https://www.seamlesswave.com/LISFLOOD8.0).

## Contacts

Github author: Martin Nguyen martinnguyen20494@gmail.com



