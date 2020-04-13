## Welcome to the companion repository for the manuscript entitled *Identifying substitution saturation in phylogenomic data*
### by David A. Duchene, Niklas Mather, Cara Van Der Wal, Simon Y.W. Ho

This repository contains the scripts used to make simulations (requires at least the R packages *ape* and *phangorn*). The main files are *satsim.R* and *satsimcont.R*. These functions simulate data using options of number of taxa, number of alignment sites, total tree length, phylogenetic tree imbalance, stemminess, and simulation model, as described in the manuscript.

The functions simulate an alignment with the desired qualities and run IQ-TREE using the desired substitution model. They return the true tree used as input (which in the study was taken from the *trees.Rdata* file), the estimated tree, and the estimated *t* entropy statistic, among other information about the run.

*satsim.R* simulates data from discrete input values, while *satsimcont.R* simulates data from a random set of values. Random values are drawn from uniform distributions for each of tree length, stemminess, and sequence length, and drawn from a set for number of taxa (8, 16, or 32), phylogenetic tree imbalance (fully balanced versus fully imbalanced), model used for simulation (JC or GTR+G).
