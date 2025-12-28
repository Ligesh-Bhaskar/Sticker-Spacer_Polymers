# Sticker-spacer polymers in active bath
This repository contains LAMMPS input scripts, data files, and generated outputs for simulations of sticker-pacer polymers immersed in passive or active Brownian particle baths. The analysis script and processed data used for figure generation are provided in the `PostProcessing/` directory, organized into `Cluster_Data/`, `Partition_Data/`, `Rg_Data/`, and `Plots/`.

### Simulations

This repository contains LAMMPS input scripts, data files, and generated outputs used to simulate sticker–spacer polymers immersed in passive or active Brownian particle baths.

- `Polymer_Stickers_ABP1200_Hybrid.dat`: LAMMPS data file defining the initial configuration of the sticker–spacer polymers and Brownian particle bath, including particle types, bonds, and simulation box dimensions.

- `Poly_Relax_ABP0.1_T300_E6_Pe0.in` \& `Poly_Relax_ABP0.1_T300_E6_Pe1.0.in`: LAMMPS input scripts used to relax the system with a passive (Pe = 0) and an active Brownian particle bath (Pe = 1.0). Polymer–polymer attractions are gradually introduced to obtain an equilibrated initial configuration.

- `Poly_Record_ABP0.1_T300_E6_Pe0.in` \& `Poly_Record_ABP0.1_T300_E6_Pe1.0.in`: Production run scripts for simulations with a passive (Pe = 0) and an active Brownian particle bath (Pe = 1.0). The outputs generated are processed to compute the average droplet size, partition coefficient, and radius of gyration.

- `Out_Relax/`: Directory containing relaxation trajectories and restart files generated during equilibration runs.

- `Out_Record/`: Directory containing production output files, and restart files used for post-processing and analysis.


### PostProcessing

The `PostProcessing/` directory contains the post-processing and plotting scripts, along with the processed data, used to generate the figures in the article:

- `Final_Plotting.py`: Script used to post-process simulation data and generate the plots used in the paper. It processes data from the directories listed below and saves the resulting figures in the `Plots/` directory.

- `Cluster_Data/`: Droplet (cluster) size as a function of Pe. This data is used to generate Fig. 5B.

- `Partition_Data/`: Partition coefficient data for different Pe values. This data is used to compute the average partition coefficient shown in Fig. 5C.

- `Rg_Data/`: Radius of gyration of the polymer in the dense and dilute phases for lower and higher binding affinities ($\epsilon = 4,\, 6$) at Pe = 0 and Pe = 0.5. This data is used to generate Fig. S5.

- `Plots/`: Final plots used in the paper.


