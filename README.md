# EHL
MATLAB solver for Elastohydrodynamic Lubrication (EHL) problems. The solver is based on a Finite Volume (FV) discretization of the Reynolds equation incorporating mass-conserving cavitation through the cavity fraction and elastic deformation through the application of the Boundary Element Method (BEM) to an elastic half-space. Shear thinning, Roelands and Dowson-Higginson relations are also embedded in the code. Furthermore, the load-balance equation is considered.

The code is split up into a steady and and unsteady solver. While the steady solver is already finished, the unsteady solver is still under construction.

The code consists out of three scripts which must be executed consecutively:

The scripts denoted with EHL_01_setup will set up the input information for the EHL solver.

The script denoted with EHL_02_mainprocess incorporates the actual EHL solver and computes the solution of the EHL problem.

The script denoted with EHL_03_visualization displays plots of the results obtained by the solver.

Furthermore, there are three reference cases:

The codes denoted with Krupka replicates the setup of:
Křupka, I., Hartl, M., Urbanec, L., & Čermák, J. (2007). Single dent within elastohydrodynamic contact-comparison between experimental and numerical results. Proceedings of the Institution of Mechanical Engineers, Part J: Journal of Engineering Tribology, 221(6), 635-644.

The codes denoted with Mourier replicates the setup of:
Mourier, L., Mazuyer, D., Lubrecht, A. A., & Donnet, C. (2006). Transient increase of film thickness in micro-textured EHL contacts. Tribology International, 39(12), 1745-1756.

The codes denoted with Rosenkranz replicates the setup of:
Rosenkranz, A., Szurdak, A., Gachot, C., Hirt, G., & Muecklich, F. (2016). Friction reduction under mixed and full film EHL induced by hot micro-coined surface patterns. Tribology International, 95, 290-297.

This code is free to use by anyone. Still, citing this repository is greatly appreciated if it is of use to your work.

Erik Hansen
