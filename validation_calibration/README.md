These folders are used to validate and calibrate the LISFLOOD-FP flood model. Download all these three folders into your local computer and start to run the scripts from the ```script``` folder. The scripts were written in jupyter notebook, so we suggest to use jupyter notebook with the kernel of the environment installed earlier. Please change the path in ```basepath```, ```os.chir```, and ```main_dir``` to the ```data``` folder path as it stores all the data sources. After that, please run the script as the order below. **_Notice_**: Land (51559) and bathymetry (50554) files were downloaded on 15/04/2024. The OpenStreetMap data was downloaded into cache folder on 25/01/2023. Hence, the reproducible work could be different because the versions of these data could be changed.

1. Run ```DEMforvalidation_instruction.ipynb```, it will create a folder named ```dem_generation``` to store all the DEM data.
2. Run ```DEMforvalidation_generation.ipynb```, it will create the DEM and all data used to create this DEM.
3. Run ```Roughnessforvalidation.ipynb```, it will create the roughness length in this folder. After this, copy the DEM and roughness length files outside of ```dem_generation``` folder to ```data``` folder.
4. Run ```Rainfall_resampling.ipynb```, it will resample the rain on grid to higher resolution.
5. Run ```Rainfall_netcdf_generation.ipynb```, it will create format file of rainfall that should be used for the LISFLOOD-FP flood model.
6. Run ```Validation_otherfiles.ipynb```, it will create necessary files to validate the data.
7. Run ```Calibration_otherfiles.ipynb```, it will create necessary files to calibrate the data.
8. Run ```Evaluationfor_validation_calibration.ipynb```, it will create the plots to visualise validation and calibration results.

The results were shown as below

https://github.com/Martin20494/Grid_Orientation/blob/main/validation_calibration/results/validation_calibration.jpg


