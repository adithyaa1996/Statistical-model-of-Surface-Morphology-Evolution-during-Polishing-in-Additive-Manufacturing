# Statistical model of Surface Morphology Evolution during Polishing in Additive Manufacturing

This GitHub repository maintains data associated with our accepted paper in IISE Transactions titled "Statistical and Dynamical model of Surface Morphology Evolution during Polishing in Additive Manufacturing". To briefly summarize,

**1. Polishing_stagewise_data.zip** - Contains height values measured at 32 different locations on the 3D printed sample using an optical profilometer prior to polishing (Stage 0) and post every stage of polishing (Stages 1 to 6). Please refer to the following paper for experimentation details and process parameters: "_Jin, S., A. Iquebal, S. Bukkapatnam, A. Gaynor, and Y. Ding (2019, 10). A gaussian process model-guided surface polishing process in additive manufacturing. Journal of Manufacturing Science and Engineering 142, 1–17._"

**2. Initial_surface_generation.m** - Script containing the Initial surface generation algorithm using the random circle packing algorithm. This file generates the surface asperity distribution and their graph connectivity of a 3D printed sample prior to polishing. One such realization is stored and compared with experimental data (Refer #3).

**3. Stage0_fitted_data.mat** - .mat file containing data pertaining to height measures of the 3D printed sample prior to polishing and generated initial surface (simulation) which is statistically similar to the actual data. 

**4. Parameter_fitting_Polishing.m** - Script containing the model capturing polishing dynamics with network formation, evaluated at each stage of polishing. This file generates the Bearing Area Curves of the initial surface generated after each stage of polishing. (It uses other functions defined in #5).

**5. surface_roughness.m, graph_evolution.m, solve_for_d.m, KLDiv.m** and **Gen_hurst.m** - Matlab scripts containing functions that are called within the main script (Parameter_fitting_Polishing.m)

**6. Simulated_Annealing.zip** - Zip file containing files related to Simulated Annealing Algorithm. Please read the **README_Simulated_Annealing.txt** for instructions to reproduce the optimized parameter solutions. 

**7. pub_fig.m** - Script containing the formatting options for plots and figures. 
