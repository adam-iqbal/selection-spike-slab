# Spike-and-slab priors for variable selection in sample selection problems
Basic implementation of the methods described in "Bayesian variable selection in sample selection models using spike-and-slab priors".

The folder "methods" includes the spike-and-slab sampler implementation, alongside implementations for stepwise selection and Adaptive LASSO. The supported priors for the spike-and-slab sampler are normal, Laplace, and t-distribution with 3 degrees of freedom, though this can be modified in the code itself.

Required packages: sampleSelection, MASS, rtruncnorm (for spike-and-slab only), numDeriv (for Adaptive LASSO only), ssmrob (for the ambulatory data only)

The folder "testing" contain parallelized functions to store results from simulations, but this needs to be run manually. An example is provided in gibbs_simulation.R for one scenario. It also contains code for the two real data applications.
In ambulatory_testing.R and rand_testing.R, the provided code will store all the output from the four models, provided all the packages and required functions are loaded into the workspace.
"speed_test.R" contains the code used to test the running time of individual iterations (which will differ by computer, of course).

An example of a simulation run can be found in "testing\gibbs_simulation.R". All simulation studies were run with the hyperparameter elicitation as in the example, with only n, p, corr and the given distribution (alpha_spike, alpha_slab, beta_spike, beta_slab) changing between runs, so choosing corresponding values for these parameters should reproduce the exact results in the paper.
