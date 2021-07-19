# Introduction
This code implements the L2-NSGA algorithms used to optimize the refractive index distribution in mode de-multiplexging devices for Hermite-Gauss beams. 
Reference papers:
 
[1] G. Di Domenico, D. Weisman, A. Panichella, D. Roitman and A. Arie, "Planar on-Chip Mode Sorter" (2021).

[2] M. J. G. Olsthoorn, A. Panichella, "Multi-objective Test Case Selection Through Linkage Learning-based Crossover",  13th Symposium on Search-based Software Engineering (2021).

# How to use the code
The code is written in `Matlab` and tested with Matlab version V9.9.0.1495850 (R2020b) Update 1. To use the code, you simply need to:
1. Open the folder in Matlab and set it as the working directory
2. Either run the file `main2.m` or `main3.m` for the two-objective and three-objective formulation of the problem, respectively.

The folder contains the following files:

* `main3.m`: main script to run the two-objective version of the problem
* `main2.m`: main script to run the three-objective version of the problem
* `BlackBoxfunctionTwoBeams.m`: objective functions to optimize with L2-NSGA for the two-objective problem
* `BlackBoxfunctionThreeBeams.m`: objective functions to optimize with L2-NSGA for the three-objective problem
* `evaluation.m`: wrapper used to evaluate an entire population (for L2-NSGA) at once
* `dominance.m`: the function that computes the non-dominance among generated solutions		
* `L2NSGA.m`, the main search algorithm which implements a linkage-learning variant of NSGA-II  as proposed by M. J. G. Olsthoorn and A. Panichella. This variant generates offspring by using the linkage structures are derived using the agglomerative hierarchical clustering.
* `extractFOS.m` function that extracts the linkage structures using the agglomerative hierarchical clustering
* `tournament_selection.m`, `non_domination_sort_mod.m`, and `replace_chromosome.m` implement the key routines of NSGA-II. These three functions are implemented by Aravind Seshadri and public available in [MathWork](https://www.mathworks.com/matlabcentral/fileexchange/10429-nsga-ii-a-multi-objective-optimization-algorithm).