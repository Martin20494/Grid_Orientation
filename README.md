# Grid_Orientation
Quantifying uncertainty in flood predictions caused by transforming square grid orientation in LiDAR-derived DEM process

## Summary

Digital elevation datas, or DEMs, are critical for generating reliable flood predictions. The most common method to generate DEMs involves sampling and interpolating LiDAR data onto a North-South square grid. However, the orientation of this grid, which can introduce variability in elevation data and thus influence flood predictions, is frequently overlooked. In this study, we quantify the uncertainty in flood model outputs caused by transforming this grid orientation. Using a Monte Carlo method, we produced multiple DEMs by randomly rotating and/or translating the square grid orientation to predict floods for uncertainty analysis. This Monte Carlo framework was also applied at various resolutions (2, 5, 10, and 20 meters) and flood return periods (5, 10, 20, 50, and 1000 years). 

Results show that the highest uncertainty in flood predictions occurred predominantly at the flood extent boundaries and near the river. Rotating the grid orientation caused more uncertainty than translating. The highest uncertainty was found when both rotating and translating the grid. The variability in both flooded areas and the number of flooded buildings followed these patterns. Finer resolutions revealed fewer variations in flood predictions and less variability in flooded areas and the number of flooded buildings. These variations were also determined by the amount of river discharge. Significant variations occurred depending on whether the river discharge was insufficient to cover the floodplain surfaces. Also, the flooded areas and the number of flooded buildings could expand or contract based on whether the water reached new locations and these places included residential buildings.


## Workflow

![transformation_process](https://github.com/Martin20494/Grid_Orientation/assets/55137629/4ecdd3b5-2e28-41b2-ae65-8ba4044b20d8)




