# Integrative ITR simulations

The execution unit is **one setting × one condition**. Open that condition's `run.R` to run its selected methods. There is no root script that launches all conditions.

The repository contains **Setting 1 (constant confounding)** and **Setting 2 (quadratic confounding)**, each with the baseline and seven perturbations. All 16 condition directories and their empty `output/` folders are provided. No simulation results are included.

## Quick start

Install the two R packages once:

```r
install.packages(c("Rcpp", "randomForest"))
```

Rcpp requires a working C++ compiler. On macOS, run `xcode-select --install` in a terminal if Apple's Command Line Tools are not already installed.

Edit the shared **setup.R** at the repository root:

```r
replications <- 1L
start_seed <- 1001L
cores <- 2L
methods <- c("WY", "Target", "Pooled", "Wu2020", "IP", "CC", "CC_DRGO")
output_folder <- "output"
```

Then choose a condition. From the repository root, run only Setting 1's baseline:

```bash
Rscript setting1/baseline/run.R
```

Or run it in R:

```r
source("setting1/baseline/run.R")
```

You can also open `setting1/baseline/` as your working directory and use:

```r
source("run.R")
```

For a different condition, use its own file, for example:

```bash
Rscript setting2/tau_underfit/run.R
```

Each `run.R` explicitly names its setting and condition. No other condition is executed. You do not need to change the working directory in the script or source shared method files yourself.

## Shared execution settings

| Parameter in setup.R | Meaning |
|---|---|
| `replications` | Number of consecutive seeds for the chosen condition |
| `start_seed` | First seed; also the only seed used for the decision-boundary figure |
| `cores` | Parallel optimizer starts within a fit; methods and seeds run sequentially |
| `methods` | Methods to execute within the chosen condition; all seven by default |
| `output_folder` | Folder name within each condition, `output` by default |

The defaults give **7 fits and 1 figure per execution**, using seed **1001**. Set `replications <- 50L` to use seeds **1001:1050**, matching the paper's cluster study. That gives **350 fits and still 1 figure per condition**. Running all 16 condition files gives 112 fits at one replication or 5,600 fits at 50 replications.

The figure contains the true boundary, the observed target contrast, and the selected methods' fitted boundaries. It shows the start seed and each method's regret for that seed. Every seed contributes to the regret tables, but later seeds are not plotted.

macOS/Linux support multiple optimizer-start workers; Windows uses one. For interactive RStudio execution, set `cores <- 1L`; use terminal `Rscript` for multicore execution. `output_folder` must be a simple folder name, not an absolute path.

## Directory structure

```text
integrative_itr/
├── README.md
├── setup.R
├── methods/                     # Shared implementation; no results here
│   ├── run_case.R               # Load setup and execute one condition
│   ├── design.R                 # True functions, bases, data generation
│   ├── methods.R                # Seven methods and nuisance estimation
│   ├── fitting.R                # Connect scores to the common optimizer
│   ├── optimizer.R              # Newton/DC, initialization, convergence
│   ├── fast_backend.R           # Load the compiled functions
│   ├── fast_kernels.cpp         # Fast loss and gradient calculations
│   ├── plots.R                  # First-seed decision boundaries
│   └── simulation.R             # Seed/method loop, checkpoints, tables
├── setting1/                    # Constant confounding
│   ├── baseline/
│   │   ├── run.R
│   │   └── output/              # Empty until you run this condition
│   ├── target_large/
│   ├── trial_large/
│   ├── no_shift/
│   ├── strong_shift/
│   ├── tau_underfit/
│   ├── c_full/
│   └── policy_full/
└── setting2/                    # Quadratic confounding
    ├── baseline/
    ├── target_large/
    ├── trial_large/
    ├── no_shift/
    ├── strong_shift/
    ├── tau_underfit/
    ├── c_full/
    └── policy_full/
```

**Every condition directory contains its own `run.R` and empty `output/`.** The invisible `.gitkeep` file preserves an empty output directory in Git; it contains no results. The old `Codes/` directory, root `run.R`, and `papersetting0`–`papersetting15` directories are not used. Extract this ZIP into a fresh directory to avoid mixing the old and new layouts.

A condition's `run.R` calls `methods/run_case.R`, which reads the root `setup.R` and loads the shared code. `simulation.R` obtains data from `design.R`, calls `fitting.R`, and writes results. `fitting.R` uses `methods.R` for scores and `optimizer.R` for fitting; the optimizer uses `fast_kernels.cpp` through `fast_backend.R`. The seven statistical methods are implemented in one file, `methods/methods.R`.

## Results for one condition

After running `setting1/baseline/run.R`:

```text
setting1/baseline/output/
├── decision_boundary.png       # Start seed only, all selected methods
├── regret.csv                  # Every selected method and seed
├── summary.csv                 # Mean regret and replication Monte Carlo SE
├── paired_comparisons.csv      # Same-seed method differences
├── settings.csv                # Sample sizes, shift, basis dimensions
├── timing.csv                  # Backend, data preparation, fit/evaluation, plot times
├── session.txt                 # Configuration and R session information
├── run_config.rds              # Configuration and source/runtime fingerprint
├── status.txt
├── completion_report.txt
└── checkpoint.rds              # Internal resume data; excluded from Git
```

There are no extra method, optimizer, basis-mode, or seed subdirectories. Each condition owns its results independently. From the repository root, read them with:

```r
read.csv("setting1/baseline/output/summary.csv")
read.csv("setting1/baseline/output/paired_comparisons.csv")
```

Regret is evaluated on 100,000 independent target covariates:

$$
\widehat{\mathcal R}(\widehat d)
=\frac{1}{N}\sum_{i=1}^{N}|\tau(X_i)|
\mathbf{1}\{\widehat d(X_i)\ne\mathbf{1}[\tau(X_i)>0]\}.
$$

Lower regret is better. `test_MCSE` measures test-sample integration error for one fitted policy. `replication_MCSE` measures variation across replications; it is `NA` with one replication. In paired tables, a negative `mean_A_minus_B` favors method A. Available numerical fits, including `optimizer_warning` fits, enter summaries; their counts are reported explicitly.

The checkpoint stores all per-seed regret records and only the start seed's policy coefficients. Later policies are discarded after evaluation. CSVs, figures, and session information can be committed with the code; `.gitignore` excludes checkpoints, compiled files, and R session clutter.

## Settings and perturbations

Let $B_q(x)$ contain normalized products of Legendre polynomials with total degree at most $q$, including the intercept. Its dimension is $(q+1)(q+2)/2$. Write $d_D$, $d_\tau$, and $d_c$ for the policy, working treatment-effect, and working confounding basis dimensions.

Both settings use $r^2=x_1^2+x_2^2$ and

$$
\tau(x)=K(r^2-0.2)(r^2-0.65)(r^2-1.2),
$$

where $K$ is the reciprocal IQR of the unscaled effect on a fixed $601\times601$ grid over $[-1,1]^2$. The true decision boundary is given by three circles, restricted to this square.

| Setting | True confounding function $c(x)$ | Baseline $(d_D,d_\tau,d_c)$ |
|---|---|---|
| Setting 1: constant confounding (`setting1`) | $2$ | $(28,28,1)$ |
| Setting 2: quadratic confounding (`setting2`) | $2+x_1-x_2+x_1^2-x_2^2$ | $(28,28,6)$ |

Target covariates are uniform on $[-1,1]^2$. Trial covariates have density proportional to $\exp[-a(x_1^2+x_2^2)]$ on the same square. Let $S=0$ denote target data and $S=1$ trial data. Independently of covariates, $U$ is equally likely to be $-1$ or $1$. Target treatment is $A=\mathbf{1}\{U=1\}$; trial treatment is independent Bernoulli$(1/2)$. The outcome is

$$
Y=x_1+x_2+(A-1/2)\tau(X)+(1-S)Uc(X)/2+\varepsilon,
\qquad \varepsilon\sim N(0,(1/3)^2).
$$

Thus the observed target contrast is $\tau+c$, whereas the trial contrast is $\tau$.

| Condition | Change from baseline, applied to both settings |
|---|---|
| `baseline` | $n_0=8000$, $n_1=800$, $a=0.75$; dimensions above |
| `target_large` | $n_0=16000$ |
| `trial_large` | $n_1=1600$ |
| `no_shift` | $a=0$ |
| `strong_shift` | $a=2$ |
| `tau_underfit` | Treatment-effect basis $B_2$: $d_\tau=6$ |
| `c_full` | Confounding basis $B_6$: $d_c=28$ |
| `policy_full` | Policy basis $B_7$: $d_D=36$ |

Only the listed component changes. True functions remain fixed across perturbations. Marginal-outcome and arm-outcome nuisance regressions are **estimated using $B_6$, dimension 28, throughout**. They are not given the oracle linear basis for $\mu_s=x_1+x_2$. In `tau_underfit`, the working effect basis alone shrinks to dimension 6; nuisance bases remain dimension 28. In `policy_full`, only the policy basis expands to dimension 36. This matches the completed cluster study.

## Methods and numerical implementation

| Code | Role |
|---|---|
| `WY` | Transport-weighted trial augmented outcome score |
| `Target` | Observational target-only augmented outcome score |
| `Pooled` | Pooled trial/target augmented outcome score |
| `Wu2020` | Target-learned random-forest features with trial outcome calibration |
| `IP` | Integrated treatment-effect plug-in policy |
| `CC` | Target score corrected by estimated pointwise confounding |
| `CC_DRGO` | Target score corrected by the DRGO estimating equation |

All methods use the common policy basis for a given case. Five-fold cross-fitting, nested nuisance fitting where required, and Wu2020's 500-tree forests follow the cluster implementation. The marginal outcome is estimated from outcomes, not supplied as its true function.

The compiled loss/gradient calculations, Newton steps, line search, and retry rules are preserved. Each fit uses **11 starts**: a zero start, nine random starts, and a scaled least-squares start. Selection uses the largest empirical policy value among those starts. The base limits are 1,000 outer DC iterations and 200 inner Newton iterations, with BFGS fallback and extended retries for unresolved starts. The `max_dc` argument in score helpers is retained for compatibility; the active limits come from `optimizer_control()` in `methods/optimizer.R`.

`completed` requires convergence of all starts for every requested fit and successful plots. `mixed` indicates remaining fit/optimizer/plot issues; `failed` means no valid numerical fits. Inspect `regret.csv` for errors and convergence diagnostics. Stationarity does not establish global optimality. True effects are used for evaluation and an after-fitting diagnostic, never for choosing starts or fitting the policy.

## Runtime and speed

The directory structure does not change the statistical calculation. The execution code avoids repeated work:

- The 100,000-point test basis is computed once per seed and shared by all methods.
- Policy normalization and the scaled design matrix are computed once per fit and shared by all 11 starts and their retries.
- Newton iterations reuse an already computed margin when evaluating the gradient.
- The routine runner disables the in-memory exact-problem cache: each method's score is solved once, so serializing the full matrix to generate an unused cache key only adds overhead. Completed fits are still reused from checkpoints.
- Compiled C++ code is cached across R sessions outside the repository, under R's user cache directory for `integrative_itr`. The cache signature includes source, R/Rcpp versions, platform, build flags, and available Makevars/Makeconf files. The first build still requires compilation. Set the environment variable `ITR_CACHE_DIR` to use a different cache location. If you replace the compiler toolchain, use a fresh cache directory.

Progress messages identify the active condition, seed, method, and plotting stage. `timing.csv` separates backend loading, data/test-basis preparation, fitting plus evaluation, and plotting. Per-method `regret.csv` includes `fit_seconds`, `evaluation_seconds`, `retried_starts`, and `total_dc_iterations` (including initial attempts). Checkpoint reuse is announced and does not rerun a fit; its stored fit times describe the original calculation, while `timing.csv` describes the current invocation.

Comparisons with older `github_full` runs need matching workloads. For example, old P0 used 4,000 target observations, 400 trial observations, and a 10-dimensional policy; the current baseline uses 8,000, 800, and 28 dimensions. Old test samples had 50,000 points; the paper uses 100,000. The supplied `github_full` optimizer also lacked the extended retries present in the paper implementation. Its `C_SCALE=2` and default seed 701 differ from this package's paper settings and seed 1001. These differences can change both runtime and results.

This revision keeps all 11 starts, convergence tolerances, retries, random-forest sizes, sample sizes, and estimators. No full-paper speedup is claimed without a matching-workload benchmark.

## Reproducibility and resuming

For matching archived cluster results, `setting1` corresponds to the historical label `S3`, and `setting2` to `S4`. This is a label-only change; functions, dimensions, conditions, seeds, and estimators are unchanged. The paper and new output directories use the consecutive numbering.

The default seed labels are the cluster study's replication indices. The original component-specific RNG offsets, fold seeds, random-forest seeds, basis column order, and optimizer starts are retained. Methods share generated data and test covariates within a case/seed. Separate RNG streams preserve paired comparisons when a sample size changes. Internally, the data object stores `seed + 1` and fitting adds its original offsets; do not add these offsets to `start_seed` yourself.

From the repository root, run one condition with single-threaded numerical libraries:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript setting1/baseline/run.R
```

R version, package versions, platform, numerical libraries, source hashes, and controls are recorded. Identical seeds do not guarantee bitwise-identical results across R/compiler/BLAS platforms. Use the recorded cluster environment when exact numerical comparisons are needed.

Re-running a condition with the same configuration reuses its completed fits. Raising `replications` from 1 to 50 adds the remaining seeds. Failed or unresolved fits are retried. Changing the start seed redraws its boundaries and, when necessary, refits that seed to recover its policy coefficients. Existing regret records are retained in checkpoints; exported tables contain only that condition, the currently selected methods, and the requested seed range. Each condition has its own independent checkpoint.

A changed code/runtime/model fingerprint causes an explicit stop instead of mixing results. Choose a different `output_folder` in `setup.R` (for example, `output_v2`) or move that condition's old output directory before starting that experiment. Optimizer-start parallelism is part of this fingerprint.

No paper simulation replications were run while preparing this revision. Validation includes R syntax checks and fixed-input numerical checks of the optimizer and evaluation path. In the checked cases, the optimized solver reproduces the previous version's coefficients, selected start, full convergence traces, and RNG state exactly, including retry and parallel-start cases. These checks do not establish an end-to-end runtime improvement for every paper condition.
