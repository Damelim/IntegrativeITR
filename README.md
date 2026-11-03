# Target-population treatment rules: simulation replication

Four files; no private source tree or downloaded data are needed. Run commands
from this directory. The paper's three 25-replication designs are S1, S2 and
S3; S2 corresponds to the archived name `S2a_quadratic_c`. S4 is a new,
paired sensitivity experiment, **not a reported result**.

## Start locally

Use R (the cluster used 4.5.3). Base R suffices for IP, CC, WY, Target,
Pooled and Fusion. The full nine-method comparison additionally needs:

~~~r
install.packages(c("randomForest", "ncvreg"))
~~~

Generate a dataset and fit a policy in R:

~~~r
source("generate_data.R")
source("methods.R")
s <- generate_data("S4", "baseline", rep = 1)
head(s$target)                         # x1, x2, A, Y; no U
head(s$trial)
cc <- fit_itr(s, "CC")
ip <- fit_itr(s, "IP")                 # exactly the same data as CC
cc$theta                               # treat if theta[1] + theta[2]*x1 > 0
c(CC = population_regret(cc$theta,"S4"),
  IP = population_regret(ip$theta,"S4"))
~~~

Or, from a terminal:

~~~bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript run.R S4 baseline CC 1        # one complete fit, not a cheaper solver
Rscript run.R S2 double_trial IP 25   # reproduce an original paper cell
Rscript run.R S4 baseline all 25      # all nine methods, one condition
~~~

`results/` contains one CSV per setting, condition, method and replication
count. To avoid overwriting previous outputs, choose a new directory with
`ITR_RESULTS=/path/to/new/results` when repeating a run. Never combine a
one-replication example with the full run: replication 1 is the same seed.

## Read the code in this order

1. **generate_data.R**: `paper_setting()` defines the paper's functions;
   `conditions` defines the design; `generate_data()` returns target, trial
   and independent evaluation samples. `population_regret()` is for evaluation only.
2. **methods.R**: `fit_itr(sample, method)` dispatches to the original fitted
   scores and DC/BFGS policy optimizer. Below that public entry point are
   the fitting kernels themselves, not imports from another repository.
3. **run.R**: `run_cell()` generates data, fits and evaluates the requested
   method. `collect_s4()` checks every cell before summarizing the experiment.

## What changes in S4?

Write x = X1, y = X2 and H = 1 + 0.35x + 0.15x² + 0.25y².

| Setting | Treatment effect τ | Confounding c | (policy, effect, correction) dimensions | Target noise SD | RCT noise SD |
|---|---|---|---|---|---|
| S1 | (x−0.25)H | 0.40x+0.10y | (2,10,2) | 1 | 1.25 |
| S2 | x(1+0.10x)H | −1.65x² | (2,15,1) | 0.04+0.65x² | 1.25 |
| S3 | (x−0.25)(1+0.20y) | 0.35+0.25x−0.25y | (2,6,3) | √(0.49+2x²) | 0.25 |
| S4 | **same as S2** | **same as S2** | **(2,15,1)** | **same as S2** | **0.625** |

In S2/S4, the optimal boundary is x=0. If the estimated correction coefficient
is in error by h, its contribution to CC's score error is hx²: at x=0 it
vanishes, and nearby it is small. IP fits the whole 15-dimensional effect.
S4 asks whether this local feature produces a **mean** advantage for CC
when trial outcomes are more precise. Halving trial noise quarters conditional
outcome variance; at a fixed RCT design it quadruples Gaussian regression
information. It benefits both estimators and does not guarantee the ranking.

Both covariates are independent zero-mean normals truncated to [−3,3].
The target SD before truncation is 1. All settings use the six conditions below.
No effect normalization, measured W, policy approximation, or known cutoff is
used in fitting. A latent fair sign U gives P(A=1|X,U,S=s)=1/2+(1−s)U/4;
Y(a)=(a−1/2)τ(X)+(1−s)Uc(X)+η_s(X)ξ_a. Thus both populations share the causal
effect, while their observed treatment contrasts are τ+c (target) and τ (RCT).

| Condition | n0 | n1 | RCT covariate SD before truncation |
|---|---:|---:|---:|
| baseline | 4000 | 100 | 1 |
| half_target | 2000 | 100 | 1 |
| double_target | 8000 | 100 | 1 |
| double_trial | 4000 | 200 | 1 |
| contraction | 4000 | 100 | 0.8 |
| expansion | 4000 | 100 | 1.25 |

## Full S4 run and output

One task = one of six conditions × one of nine methods (WY, Target, Pooled,
Wu2020, IP, CC, Fusion, CC_SCAD, Fusion_SCAD) × all 25 matched replications.
DC uses 16 starts per optimization, up to 200 DC iterations per start. The
threshold class also contains both constant policies. Do not drop a
replication because the fitted policy treats everyone or no one.

To run the equivalent of the cluster batch locally (potentially slow), in a
new results folder:
~~~bash
export ITR_RESULTS=results-full
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
for task in $(seq 1 54); do Rscript run.R task "$task" || break; done
Rscript run.R collect
~~~

The collector checks all 1,350 fits; it writes `s4_r25_replications.csv`,
`s4_r25_all_methods_summary.csv` and `s4_r25_paired_summary.csv`.
Primary comparison: baseline mean of regret(CC) minus regret(IP) with its
paired Monte Carlo standard error. Negative favors CC. Report all five
sensitivity conditions, wins, ties, median and largest loss, regardless of
direction. Weight RMSE is evaluated only for unpenalized Fusion. Population
regret integrates disagreements; the second regret column uses an independent
50,000-observation test sample.

S4 reuses S2's 25 streams (seeds 9730001–9730025) on purpose: covariates,
treatments, target outcomes and standardized RCT noise match. It is a
controlled sensitivity analysis informed by S2, **not independent validation**.
After collection, compare the change in CC-minus-IP within the same seeds:
~~~bash
Rscript run.R compare "/path/to/dc9_main3_r25_replications.csv"
~~~
Cross-platform R/BLAS/package differences can change numerical optimizer
paths. Use the original cluster CSVs as the record of the paper results;
this release does not overwrite them. To audit one archived baseline fit
without writing any results, use:
~~~bash
Rscript run.R verify "/path/to/dc9_main3_r25_replications.csv" S2 CC 1
~~~
This prints both values and their differences, including the fitted policy
coefficients and optimizer call counts. For a strict comparison, run on the
same R installation, package versions and BLAS as the original.
