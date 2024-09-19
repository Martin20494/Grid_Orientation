# Grid Alignment

### Title

Quantifying uncertainty in flood predictions LISFLOOD-FP due to arbitrary conventions in grid alignment in LiDAR-derived digital elevation models

### Authors

Martin Nguyen, Prof. Matthew D. Wilson, Dr. Emily M. Lane, Prof. James Brasington, and Dr. Rose A. Pearson

## Overview

Digital elevation models (DEMs) generated by sampling and interpolating LiDAR data onto a square grid can produce reliable flood predictions. However, the arbitrary conventions in grid alignment that can introduce uncertainty in flood predictions are frequently overlooked. Hence, our research quantified this uncertainty using a Monte Carlo approach and flood model LISFLOOD-FP to generate multiple flood simulations for analysis. The results demonstrate 7% ($\approx$ 168000/2384476) variability in flood extent and 27% ($\approx$ 50/185.94) variability in the number of flooded impacted by this uncertainty.

<div align="center">
	<img width = "45%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data/All_results/Results_representative/S3_area.jpg">
	<img width = "45%" src="https://github.com/Martin20494/Grid_Orientation/blob/main/other_files/data/All_results/Results_representative/S3_building.jpg">
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



