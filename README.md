# Sticker-spacer polymers in active bath
This repository contains LAMMPS input scripts, data files, and generated outputs for simulations of sticker-pacer polymers immersed in passive or active Brownian particle baths. The analysis script and processed data used for figure generation are provided in the Analysis directory, organized into Cluster_Data, Rg_Data, Partition_Data, and Plots.

### Analysis

The `Analysis/` directory contains the post-processing and plotting scripts, along with the processed data, used to generate the figures in the article:

- `Final_Plotting.py`: Script used to post-process simulation data and generate the plots used in the paper. It processes data from the directories listed below and saves the resulting figures in the `Plots/` directory.

- `Cluster_Data/`: Droplet (cluster) size as a function of Pe. This data is used to generate Fig. 5B.

- `Partition_Data/`: Partition coefficient data for different Pe values. This data is used to compute the average partition coefficient shown in Fig. 5C.

- `Rg_Data/`: Radius of gyration of the polymer in the dense and dilute phases for lower and higher binding affinities ($\epsilon = 4,\, 6$) at Pe = 0 and Pe = 0.5. This data is used to generate Fig. S5.

- `Plots/`: Final plots used in the paper.


