# mtdna-network-partition
Simulation and statistics of mtDNA behaviour at cell divisions

Requirements
====
Bash for wrapper scripts, C, R with `ggplot2` and `gridExtra` for visualisation

Wrapper scripts
====
* `run.sh` -- wrapper script for simulation
* `analyse-all.sh` -- wrapper script for analysis. Run this after the jobs from `run.sh` have terminated.
* `run-rev.sh` -- wrapper script for simulation with reduced number of simulations, for review.

Code
====
* `network-sim.c` -- code to simulate partitioning of mtDNAs given some simulated network structure. Takes nine command-line parameters:

  * switches: --snapshot or --simulate -- whether to output cell snapshot only, or to simulate the full number of simulations. 
  * prop: the proportions of the parent cell apportioned to the smallest daughter -- ranges from 0 to 0.5 noninclusive.
  * het: initial heteroplasmy -- ranges from 0 to 1 inclusive
  * nseeds: nonnegative integer -- number of intial growth points for the network.
  * nsim: nonnegative integer -- number of simulations.
  * p and q: the proportions of wildtype mtDNA and mutant mtDNA, respectively, to include in the network -- both range from 0 to 1 inclusive.
  * lambda: ranges from 0 to 0.1 inclusive -- whether to perturb mtDNA molecules before outputing statistics.
  * repel: 1 (true) or 0 (false) -- whether mtDNA molecules are self-avoiding or not

For example `./network-sim.ce --snapshot 0.5 0.1 4 0.5 0.5 0.04 0` outputs summary statistics and simulation snapshots for het = 0.1, nseeds = 4, p=q=0.5, lambda = 0.04 and no self-avoidance of mtDNA molecules.

* `plotcell.R` -- plots visualisations from simulation snapshots
* `compare-both-new.R` -- compares physical simulation, statistical simulation, and theory for different statistics
* `plotlambda.R` -- plots influence of post-positioning diffusion
