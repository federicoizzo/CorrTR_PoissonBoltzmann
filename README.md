# CorrTR_PoissonBoltzmann

Corrected Trapezoidal Rule for the Poisson-Boltzmann problem in 3D

Open the folder, start Julia, and run include("initialize.jl") to initialize the setting. 
There will be packages to install if not present, specifically:
- Dates
- DelimitedFiles
- SpecialFunctions
- PyPlot
- QuadGK
- LinearAlgebra
- Statistics

Then, you can open corrected_TR_main.jl to modify the scripts you wish to run. 
They use the function defined in corrected_TR_fun.jl which sets up the 3D problem and computes the potentials.
In corrected_TR_sub.jl you find all the accessory subroutines, from the layer kernels definitions to the finite difference schemese to the closest point mapping setup functions to the weight computation and approximation to the matrix-vector function corresponding to the potential evaluation.
