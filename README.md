GPhase: Greedy Approach for Accurate Haplotype Inferencing
==================

Implementation for GPhase algorithm for haplotype inferencing, reported in the paper: "GPhase: Greedy Approach for Accurate Haplotype Inferencing", Kshitij Tayal, Naveen Sivadasan, Rajgopal Srinivasan.

GPhase_Parameter.ipynb - Calculation of Model Parameters used by GPhase, namely mutation rate , quadrature points and quadrature weights.This interactive ipython notebook is just for playing around with the parameters. 

Minheap.py - Min Heap implementation for maintaining top-k largest solutions. This file is used by GPhase.py and GPhase_LOOC.py.

GPhase_LOOC.py - GPhase implementation for leave-one-out cross-validation (LOOCV) where known set of A<sub>n</sub> is created by including all haplotypes except one pair and we infer the constituent haplotype from the combined genotype of remaining pair. Usage: GPhase_LOOC.py <arg1> <arg2> . Collection of haplotypes is specified as Arg1 as a comma separated file where each row contains a haplotype, each consisting of l loci.  

GPhase.py
