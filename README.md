# EHL-FBNS Solver
MATLAB solver for unsteady Elasto-Hydrodynamic Lubrication (EHL) problems, based on the EHL-FBNS algorithm. The solver uses the Finite Volume (FV) discretization of the Reynolds equation incorporating mass-conserving cavitation through the cavity fraction and elastic deformation through the application of the Boundary Element Method (BEM) to an elastic half-space. Shear thinning, Roelands and Dowson-Higginson relations are also embedded in the code. Furthermore, the load-balance equation is considered.


The code consists out of three main scripts which must be executed consecutively and two optional scripts:  
Optional: the scripts denoted with EHL_00_run will set up and run severeal EHL solver simulations within one study.  
The scripts denoted with EHL_01_setup will set up the input information for the EHL solver.  
The script denoted with EHL_02_mainprocess incorporates the actual EHL solver and computes the solution of the EHL problem.  
The scripts denoted with EHL_03_visualization displays plots of the results obtained by the solver.  
Optional: the script denoted with EHL_04_animation creates videos of unsteady simulations.


Furthermore, there are four simulation studies:  
The codes denoted with Study_A replicate the unsteady setup of:  
Mourier, L., Mazuyer, D., Lubrecht, A. A., & Donnet, C. (2006). Transient increase of film thickness in micro-textured EHL contacts. Tribology International, 39(12), 1745-1756.  
The codes denoted with Study_B replicate the steady setup of:  
Bertocchi, L., Dini, D., Giacopini, M., Fowell, M. T. & Baldini, A. (2013). Fluid film lubrication in the presence of cavitation: a mass-conserving two-dimensional formulation for compressible, piezoviscous and non-Newtonian fluids Tribology International. Elsevier, 67, 61-71.  
The codes denoted with Study_C4 replicate the steady setup of:  
Woloszynski, T.; Podsiadlo, P. & Stachowiak, G. W. Efficient Solution to the Cavitation Problem in Hydrodynamic Lubrication Tribology Letters, 2015, 58.  
The codes denoted with Study_E4 replicate a steady setup for comparison to the analytical solution of:  
Fowell, M., Olver, A. V., Gosman, A. D., Spikes, H. A. & Pegg, I. (2007). Entrainment and Inlet Suction: Two Mechanisms of Hydrodynamic Lubrication in Textured Bearings. Journal of Tribology, 129, 336-347.  

This code is free to use by anyone. Still, citing this repository or the research below is greatly appreciated if it is of use to your work. 
The code documentation and research associated with Study_A, Study_A2,Study_A3, Study_B, Study_C4 and Study_E4 is published within the article:  
Hansen, E., Kacan, A., Frohnapfel, B. & Codrignani, A. (2022). An EHL Extension of the Unsteady FBNS Algorithm Tribology Letters. Springer, 70, 1-25.  
It can be found at:  
Website:                    https://link.springer.com/article/10.1007/s11249-022-01615-1  
DOI:                        https://doi.org/10.1007/s11249-022-01615-1  
Authors:                    Erik Hansen, Altay Kaccan, Bettina Frohnapfel, Andrea Codrignani  
Published:                  15 Jun, 2022 


24.06.2022, Erik Hansen

This research was funded by Deutsche Forschungsgemeinschaft (DFG) Project Number 438122912.
