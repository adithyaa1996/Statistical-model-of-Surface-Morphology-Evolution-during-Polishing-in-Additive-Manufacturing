The matlab files MValidate1.m, MValidate2.m,….,MValidate6.m serves as the main experiment files, where we use global optimization (simulated annealing)
 to find the best set of hyper parameters by running the dynamic model forward for a stage of polishing, and evaluating the KL distance between the obtained 
asperity distribution and the asperity distribution in the data. We choose the set of parameters that minimizes the objective function comprised of a convex combination 
of the KL distance and the surface roughness metric (SA). 

We find the optimal solutions returned by simulannealbnd(), which is the native Matlab implementation of simulated annealing, and write out 
the obtained solutions in the matlab files MValidate.m (in the variables p1, p2, etc and pkT1, pkT2 etc). 

There are helper functions, such as setModelHeightsSTGk, setModelRadiiSTGk, retrieveModelHeightsSTGk, 
retrieveModelRadiiSTGk, (k in {0,1,2,…,6}) that help in setting the initial values of the height and radii distributions for each stage of polishing, 
by running the dynamics forward for the previous stages using the best found hyperparameters so far.

The matlab files objectiveSTG0To1Ver3.m, objectiveSTG1To2.m, objectiveSTG2To3.m, objectiveSTG3To4.m, 
objectiveSTG4To5.m, objectiveSTG5To6.m returns the value of the objective function for one stage of polishing by 
simulating forward the dynamics parameterized by some realization of the hyperparameters. 
We use the function percentDivergence() to evaluate the KL distance.

For Stage0, run the first part of the file, “Initialize Heights, read heights”. Then use the section of the 
file “random search, stochastic edition” to find an initial set of parameters that minimizes the 
objectiveSTG0To1 function (use the input parameters to the function t). Use the set of hyper-parameters 
found to serve as a warm start for the Simulated annealing global optimizer, in the 
section “Ideas for Experiment Design”. Once you have the parameters, you can evaluate its 
performance in the section “Evaluate”. You can use it calculate the average objective value in 
the section “Average Objective Value”. You can find the difference in SA and KLDiv between 
using the sections “Difference in KLDiv” and “Difference in SA values” respectively. You can use the 
Plot section to create the bearing area curve comparisons. Use the values of KLDiv and SA to annotate the plots. 

The same overall structure follows for the next stages too. We first initialize heights and radii, 
and then do a random search to find a warm start solution, which we feed to the global optimizer. 
Once we get the solved hyperparameters, we use them to find the average objective value, and the 
KLdiv and SA difference for the next stage heights obtained experimentally and from data. 


