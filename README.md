# mtdna-network-partition
Simulation and statistics of mtDNA behaviour at cell divisions

Requirements
====
Bash for wrapper scripts, C, R with `ggplot2` and `gridExtra` for visualisation

Wrapper scripts
====
* `run.sh` -- wrapper script for simulation
* `analyse-all.sh` -- wrapper script for analysis. Run this after the jobs from `run.sh` have terminated.

Code
====
* `network-sim.c` -- code to simulate partitioning of mtDNAs given some simulated network structure. Takes two command-line parameters: initial heteroplasmy and number of network seed points, for example `./network-sim.ce 0.1 4`. Outputs summary statistics and simulation snapshots.
* `plotcell.R` -- plots visualisations from simulation snapshots
* `compare-both-new.R` -- compares physical simulation, statistical simulation, and theory for different statistics
* `plotlambda.R` -- plots influence of post-positioning diffusion
