# Grid_Orientation
Uncertainty in the procedure of developing Digital Elevation Models from LiDAR data for flood modelling

## Summary

Computational models of flood inundation allow rainfall, river discharge and tide levels to be related to inundation extent, and are useful tools in flood risk management. A full understanding of the errors which stem from various sources such as model inputs is vital because of the important and sometimes life-critical decisions that are made based on such models. Research of Ozdemir et al. (2013) about the uncertainties in grid size and friction parameters using LISFLOOD FP model, or the uncertainties in boundary conditions using HEC RAS model (Pappenberger et al., 2005) are two typical examples.

Accurate representations of topography are an important input for flood inundation model, and have not been widely studied. For flood hazard assessment, Airborne LiDAR point cloud data is sampled and interpolated onto a Cartesian grid (raster) to create a Digital Elevation Model (DEM) which is suitable for use in a flood model. Usually, grid alignment is not considered in the processing. However, considering orientation in sampling process may introduce variability in the resulting elevation model, leading to uncertainty that propagates through to flood model output. This may be particularly apparent for raster grid-based models, where the routing of water flow on the grid may not align with environmental features such as drainage channels.

This project investigated the variation in the outputs of a flood model using a Monte-Carlo procedure, where multiple, equally likely DEMs are derived from LiDAR by adjusting the alignment (rotation) and point of origin of the model grid, and each used to predict flood inundation. LiDAR data for the Waikanae River area, New Zealand,  were rotated (0 - 90 degrees) and translated (0 - 5 meters north and east) to produce 684 10 m resolution DEMs. The hydraulic model LISFLOOD-FP (Bates et al. 2010) was then used to generate the outputs within the Monte-Carlo framework. Model results were then transformed back to their original positions so that the variability could be analysed statistically, allowing the sensitivity of model predictions to grid alignment to be assessed.


## References

Ozdemir, H., Sampson, C. C., de Almeida, G. A. M., & Bates, P. D. (2013). Evaluating scale and roughness effects in urban flood modelling using terrestrial LIDAR data. Hydrology and Earth System Sciences, 17(10), 4015–4030. https://doi.org/10.5194/hess-17-4015-2013

Pappenberger, F., Matgen, P., Beven, K. J., Henry, J. B., Pfister, L., & Fraipont, P. (2006). Influence of uncertain boundary conditions and model structure on flood inundation predictions. Advances in Water Resources, 29(10), 1430–1449. https://doi.org/10.1016/J.ADVWATRES.2005.11.012 

Bates, P. D., Horritt, M. S., & Fewtrell, T. J. (2010). A simple inertial formulation of the shallow water equations for efficient two-dimensional flood inundation modelling. Journal of Hydrology, 387(1-2), 33-45.

