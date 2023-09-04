# Statistical model of Surface Morphology Evolution during Polishing in Additive_Manufacturing

This GitHub repository maintains data associated with our accepted paper in IISE Transactions titled "Statistical and Dynamical model of Surface Morphology Evolution during Polishing in Additive Manufacturing". The reproducibility report pdf includes instructions of using the data and codes. To briefly summarize, the repository contains 8 documents.

1. Polishing_stagewise_data.zip - Contains height values measured at 32 different locations on the 3D printed sample using an optical profilometer prior to polishing (Stage 0) and post every stage of polishing (Stages 1 to 6). Please refer to the following paper for experimentation details and process parameters: "_Jin, S., A. Iquebal, S. Bukkapatnam, A. Gaynor, and Y. Ding (2019, 10). A gaussian process model-guided surface polishing process in additive manufacturing. Journal of Manufacturing Science and Engineering 142, 1–17._"

2. Parameter_fitting_Polishing.m - Matlab .m file containing the model for capturing polishing dynamics with network formation, evaluated at each stage of polishing.

3. Stage0_fitted_data.mat - .mat file containing data pertaining to height measures of the 3D printed sample prior to polishing and simulated initial surface which is statistically similar to the actual data. 

4. surface_roughness.m, graph_evolution.m, solve_for_d.m, KLDiv.m and Gen_hurst.m - Matlab scripts containing functions that are called within the main script (Parameter_fitting_Polishing.m)
