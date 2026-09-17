# Selection on Moral Hazard in Health Insurance

Estimation programs for:

> Einav, Liran, Amy Finkelstein, Stephen Ryan, Paul Schrimpf, and Mark Cullen. "Selection on Moral Hazard in Health Insurance." *American Economic Review* 103, no. 1 (2013): 178–219.

The files in this repository generate the primary results in the paper. To actually run the estimation, the original data file (`al.csv`) is required. The data is proprietary and cannot be posted online. For further details about obtaining access to the data, please contact Liran Einav (<leinav@stanford.edu>) or Amy Finkelstein (<afink@mit.edu>). For questions about the code, please contact Paul Schrimpf (<schrimpf@mail.ubc.ca>).

> **Note:** The code has been updated and verified to run on current versions of GNU Octave (tested on Octave 9.x), in addition to the original targets of Matlab 7.10.0 (R2010a) and Octave 3.6.1.

---

## Synthetic Data

Because the original data is proprietary, **privacy-preserving synthetic datasets** are provided in `efrsc_code/postEst/csv/`:

| File | Description |
|------|-------------|
| `revisedCodeSynthetic1x.csv` / `.mat` | Synthetic dataset with the same number of individuals as the real data (~4,474) |
| `revisedCodeSynthetic10x.csv` / `.mat` | Synthetic dataset 10× larger (~44,740 individuals) |

The synthetic data are generated from the posterior distribution of the estimated model. Individual covariates have been perturbed to preserve privacy (discrete variables via randomized response; continuous variables via bounded noise), while maintaining the empirical support and distributional features of the original data. Observation fields (actual plan choices and spending) are replaced with simulated values drawn from the model. The CSV format is identical to the original `al.csv`.

To generate new synthetic data, run `postEst/createSyntheticData.m` from within the `efrsc_code/` directory:

```matlab
% From efrsc_code/postEst/
opt = struct();
opt.resultFile = 'revisedCode01';   % MCMC results file (without .mat)
opt.nSim       = 4474;              % number of synthetic individuals
opt.seed       = 1337;              % RNG seed
opt.csvFile    = 'csv/mySynthetic.csv';
opt.matFile    = 'csv/mySynthetic.mat';
createSyntheticData(opt);
```

---

## Running the Code

### Prerequisites

The majority of the code is for **Matlab** or **GNU Octave**. The code has been verified to run in:
- Matlab version 7.10.0.499 (R2010a)
- Octave version 3.6.1 and current Octave 9.x releases

Computationally intensive portions are written as C-language mex files. To compile them you will need:

1. **Matlab or Octave**
2. **gcc** (or another C compiler) with **OpenMP** support (tested with gcc 4.1.2, 4.6.3, and current releases)
3. **LAPACK** (included with Matlab and Octave; see `make.m` for details)
4. **MPFR** — compile with `--enable-thread-safe`<sup>2</sup>
5. **NLOPT**<sup>1</sup>

`make.m` attempts to download **ARMS**<sup>3</sup> automatically.

### Setup

1. Install Matlab or Octave
2. Install gcc (or another C compiler with OpenMP support)
3. *(Optional)* Install LAPACK — already bundled with Matlab/Octave
4. Install MPFR (configure with `--enable-thread-safe`)
5. Install NLOPT
6. Modify `make.m` to reflect your system paths
7. Compile the mex files by running `make.m`

### Estimation

Execute `runGibbs.m` to produce Markov Chain Monte Carlo estimates:

```matlab
runGibbs
```

This will take considerable time. Many options can be set in `runGibbs.m` to reproduce the robustness results in appendix table A8. The options are currently set to produce the baseline estimates.

### Graphs and Tables

After producing estimates with `runGibbs.m`:

1. Edit `postEst/config.m` to point to the correct result file name
2. Run `postEst/runcf.m`

Tables and figures are written to:
- `postEst/tex/tables/`
- `postEst/csv/`
- `postEst/figures/`

Not all generated tables and figures appear in the paper. See the file descriptions below for the correspondence.

---

## File Descriptions

### Estimation — M-files

| File | Description |
|------|-------------|
| `check.m` | Debugging utility: checks that spending and latent variables are consistent |
| `combineSigOPbeta.m` | Combines the parts of Σ sampled in various steps (appendix D) |
| `dropBad.m` | Deletes observations from the data structure |
| `gibbs.m` | Runs the main MCMC loop |
| `gqzero.m` | Computes Gauss–Hermite nodes and weights |
| `isOctave.m` | Returns `true` when called from Octave, `false` from Matlab |
| `iwishrnd.m` | Generates inverse-Wishart random matrices (not needed in Matlab, which has a builtin) |
| `loadData.m` | Reads data from a CSV file and stores it in a structure |
| `make.m` | Compiles mex files |
| `runGibbs.m` | Main estimation script |
| `setSeed.m` | Sets the RNG seed with a unified interface for Octave and Matlab |

### Estimation — C files

| File | Description |
|------|-------------|
| `alcoa.c` / `alcoa.h` | Shared functions: utility, choice computation, random number generators |
| `chol.c` | Cholesky decomposition (Octave workaround; not needed in Matlab) |
| `findChoicesNMH.c` | Compute choices under no moral hazard (spending = λ) |
| `findChoices.c` | Compute choices |
| `findChoicesMultMH.c` | Choices for a multiplicative moral hazard model *(not in paper)* |
| `findChoiceNmultMH.c` | Multiplicative MH choices when spending = λ *(not in paper)* |
| `findValidLatent3.c` | Finds initial latent-variable values consistent with observed choices |
| `findValidLatentMultMH.c` | Same as above for the multiplicative MH model |
| `kronEye.c` | Computes the Kronecker product X ⊗ I |
| `randdtn.c` | Draws truncated-normal random variables |
| `sampleLambdaOmega.c` | Samples λ and ω (appendix D) |
| `sampleLambdaOmegaMultMH.c` | Same for the multiplicative model |
| `sampleLamlo.c` | Samples κ (appendix D) |
| `sampleMu.c` | Samples μ_λ (appendix D) |
| `samplePsi.c` | Samples ψ (appendix D) |
| `sampleShape.c` | Samples γ₂ (appendix D) |
| `sampleSigL.c` | Samples σ_λ (appendix D) |
| `XtTimesKronSIdTimesY.c` | Computes X′(Σ ⊗ I)Y |

### Post-Estimation — M-files (`postEst/`)

| File | Description |
|------|-------------|
| `balance15.m` | Finds prices such that a given fraction of simulated observations choose plan 5 when only plans 1 and 5 are offered |
| `config.m` | Sets which MCMC result file to use for tables and figures |
| `createSyntheticData.m` | Generates privacy-preserving synthetic data from the estimated model; saves as `.mat` and `.csv` |
| `exCost.m` | Computes expected spending and cost |
| `findcfPrice.m` | Computes counterfactual prices (table 10) |
| `mnrnd.m` | Samples from a multinomial distribution (not needed in Matlab) |
| `mycorr.m` | Computes correlation |
| `mycov.m` | Computes covariance |
| `newRep.m` | Creates tables 8 and 9, and figure 2 |
| `p10plot.m` | Creates additional figures *(not in paper)* |
| `paperTablesCSV.m` | Creates CSV version of table 7 |
| `paperTables.m` | Creates LaTeX version of table 7 |
| `plotLatent.m` | Creates additional figures *(not in paper)* |
| `plotP5psiOm.m` | Creates figure 3 |
| `plotSpend.m` | Additional function for figure 3 |
| `randInitReport.m` | Plots MCMC diagnostics *(not in paper)* |
| `runcf.m` | Main script: creates all paper graphs and tables |
| `simulate.m` | Simulates a dataset from estimated parameters; supports synthetic/privacy-preserving mode |
| `slideTables2.m` | Creates figures 3 and 4 |
| `spending.m` | Computes spending |
| `welfare.m` | Creates table 10 |
| `writecsv.m` | Writes a matrix to a CSV file |

---

## Correspondence Between Output Files and Paper

| Paper | EPS / TEX file | CSV file |
|-------|----------------|----------|
| Figure 2 | `prefixSpendFlex.eps`, `prefixSpendSel.eps` | `prefixspendFlex.csv`, `prefixspendSelect.csv` |
| Figure 3 | `prefixP5latent.eps` | `prefixP5latent.csv` |
| Figure 4 | `prefixAXP5latent.eps` | `prefixP5latentAvgXnoCorr.csv` |
| Table 7  | `prefixParmP.tex` | `prefixParmPse.csv` |
| Table 8  | `prefixfit.tex` | — |
| Table 9  | `prefixds.tex` | `prefixselectMH.csv` |
| Table 10 | `prefixs00WelfareP.tex` | — |

> **Note:** Generated figures and LaTeX tables will not exactly match the formatting in the paper. The paper's figures and tables were created in Excel from the CSV files produced by these programs.

---

## Notes on Octave Compatibility

The code has been updated to run correctly on current versions of Octave. Key changes from the original release include:

- Bug fixes to `sampleMu.c` (panel covariance indexing) and `findChoicesNMH.c`
- Compatibility fixes in `postEst/balance15.m` and `postEst/findcfPrice.m` for non-finite value handling
- `setSeed.m` updated to support both old and new Octave RNG interfaces
- Addition of `postEst/createSyntheticData.m` for generating synthetic datasets

---

<sup>1</sup> Johnson, Steven G. (2010). "The NLopt nonlinear-optimization package." http://ab-initio.mit.edu/nlopt

<sup>2</sup> Hanrot, Guillaume, Vincent Lefèvre, Patrick Pélissier, Philippe Théveny, and Paul Zimmermann (2010). "The GNU MPFR Library." http://www.mpfr.org

<sup>3</sup> Gilks, Wally (1997). "Adaptive Metropolis rejection sampling (ARMS)." http://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/arms.method/arms_method.zip
