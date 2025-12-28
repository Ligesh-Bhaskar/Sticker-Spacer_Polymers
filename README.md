# Sticker-spacer polymers in active bath
This repository contains LAMMPS input scripts, data files, and generated outputs for simulations of sticker–-pacer polymers immersed in passive or active Brownian particle baths. The analysis script and processed data used for figure generation are provided in the Analysis directory, organized into Cluster_Data, Rg_Data, Partition_Data, and Plots.

### Analysis

The `Analysis/` directory contains the postprocessing/plotting script and processed data used to generate the figures in the article:
- `Cluster_Data/`: Size of the droplets as a function of Pe. This data is used for generating Fig. 5B.
- `Partition_Data/`: Partition coefficients data for different Pe. This data is used for computing the average partition coefficient shown in Fig. 5C.
- `Rg_Data/`: Radius of gyration of the polymer in dense and dilute phase for lower and higher binding affinities ($\epsilon = 4, \, 6$) for Pe = 0 and Pe =0.5. It is plotted in Fig. S5.
 
- `Plots/`: final plots used in the paper

