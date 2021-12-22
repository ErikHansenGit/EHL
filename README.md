# EHL-FBNS Solver
MATLAB solver for unsteady Elasto-Hydrodynamic Lubrication (EHL) problems, based on the EHL-FBNS algorithm. The solver uses the Finite Volume (FV) discretization of the Reynolds equation incorporating mass-conserving cavitation through the cavity fraction and elastic deformation through the application of the Boundary Element Method (BEM) to an elastic half-space. Shear thinning, Roelands and Dowson-Higginson relations are also embedded in the code. Furthermore, the load-balance equation is considered.


The code consists out of three main scripts which must be executed consecutively and two optional scripts:

Optional: the scripts denoted with EHL_00_run will set up and run severeal EHL solver simulations within one study.

The scripts denoted with EHL_01_setup will set up the input information for the EHL solver.

The script denoted with EHL_02_mainprocess incorporates the actual EHL solver and computes the solution of the EHL problem.

The scripts denoted with EHL_03_visualization displays plots of the results obtained by the solver.

Optional: the script denoted with EHL_04_animation creates videos of unsteady simulations.


Furthermore, there are two simulation studies:

The codes denoted with Study_A replicate the unsteady setup of:
Mourier, L., Mazuyer, D., Lubrecht, A. A., & Donnet, C. (2006). Transient increase of film thickness in micro-textured EHL contacts. Tribology International, 39(12), 1745-1756.

The codes denoted with Study_C4 replicate the steady setup of:
Woloszynski, T.; Podsiadlo, P. & Stachowiak, G. W. Efficient Solution to the Cavitation Problem in Hydrodynamic Lubrication Tribology Letters, 2015, 58.


This code is free to use by anyone. Still, citing this repository or the research below is greatly appreciated if it is of use to your work.

Furthermore, the code documentation and research associated with Study_A and Study_C4 is submitted within the article "An EHL Extension of the Unsteady FBNS Algorithm" to Tribology Letters and is currently under review. The preprint can be found at Research Square:

Website:                    https://www.researchsquare.com/article/rs-1064824/v1

DOI:                        https://doi.org/10.21203/rs.3.rs-1064824/v1

Authors:                    Erik Hansen, Altay Kaçan, Bettina Frohnapfel, Andrea Codrignani

Submission date:            09 Nov, 2021

Posted on Research Square:  20 Dec, 2021


22.12.2021, Erik Hansen
