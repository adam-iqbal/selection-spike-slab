# Bayesian Variable Selection in Sample Selection Models Using Spike-and-Slab Priors

## Overview

This repository contains the R code used to reproduce the results in:

> Iqbal, A., Ogundimu, E.O. and Rubio, F.J. (2026). Bayesian variable selection
> in sample selection models using spike-and-slab priors.
> *Computational Statistics*, in press.
> [DOI: 10.1007/s00180-026-01745-3](https://doi.org/10.1007/s00180-026-01745-3) |
> [Preprint (arXiv:2312.03538)](https://doi.org/10.48550/arXiv.2312.03538)

**Sample selection models** (also known as Heckman-type models) arise when the
observed outcome is subject to a selection mechanism — for example, wages are
only observed for individuals who choose to work. Standard variable selection
methods do not account for this two-equation structure. This paper proposes a
**spike-and-slab Gibbs sampler** for joint Bayesian variable selection in both
the selection and outcome equations of a sample selection model. Stepwise selection and Adaptive
LASSO are also implemented as frequentist benchmarks.

## Requirements

Install the required R packages before running any script:

```r
install.packages(c("sampleSelection", "MASS", "rtruncnorm",
                   "numDeriv", "ssmrob"))
```

| Package | Required for |
|---|---|
| `sampleSelection` | All methods |
| `MASS` | All methods |
| `rtruncnorm` | Spike-and-slab sampler only |
| `numDeriv` | Adaptive LASSO only |
| `ssmrob` | Ambulatory data application only |

## Repository structure

```
selection-spike-slab/
├── methods/
│   ├── gibbs_sampler.R        # Spike-and-slab Gibbs sampler (main method)
│   ├── stepwise.R             # Stepwise variable selection (benchmark)
│   └── adaptive_lasso.R       # Adaptive LASSO (benchmark)
└── testing/
    ├── gibbs_simulation.R     # Example simulation run (one scenario)
    ├── ambulatory_testing.R   # Real data application: ambulatory care data
    ├── rand_testing.R         # Real data application: RAND health insurance data
    └── speed_test.R           # Running time benchmarks (machine-dependent)
```

> The exact filenames within `methods/` and `testing/` should be verified
> by browsing the repository, as the names above are inferred from the README
> descriptions.

## Usage

### Spike-and-slab sampler

The main function is the Gibbs sampler in `methods/`. Three slab prior
distributions are supported: **normal**, **Laplace**, and **t(3)**. The prior
type can be set via the relevant argument; other hyperparameters
(`alpha_spike`, `alpha_slab`, `beta_spike`, `beta_slab`) and sampler settings
can be adjusted via comments in the script.

### Reproducing the simulation study

An example of a single simulation scenario is provided in
`testing/gibbs_simulation.R`. All simulation studies in the paper use the
same hyperparameter elicitation as this example, with only `n`, `p`, `corr`,
and the prior hyperparameters varying between runs. Matching these values
should reproduce the exact results reported in the paper.

> **Note:** the parallelised simulation functions in `testing/` store results
> automatically but must be run manually.

### Reproducing the real data applications

Run `testing/ambulatory_testing.R` and `testing/rand_testing.R` after loading
all required packages and functions into the workspace. Each script stores the
output from all four models.

## Related resources

- [bvsss-nlp](https://github.com/adam-iqbal/bvsss-nlp) — extension to
  Bayesian variable selection under sample selection and model misspecification
- [BVSSurv](https://github.com/FJRubio67/BVSSurv) — related short course on
  Bayesian variable selection for survival data

## Citation

If you use this code, please cite:

```bibtex
@article{iqbal2026spikeslab,
  author  = {Iqbal, A. and Ogundimu, E.O. and Rubio, F.J.},
  title   = {Bayesian variable selection in sample selection models using
             spike-and-slab priors},
  journal = {Computational Statistics},
  year    = {2026},
  note    = {in press},
  doi     = {10.1007/s00180-026-01745-3}
}
```
