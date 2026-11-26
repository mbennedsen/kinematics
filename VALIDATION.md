# Validation

The package provides 3 main figures, 10 supplementary figures and 23 supplementary tables. Validation covers the archived numerical results, reproduction of all 36 items, selected fresh estimates, file integrity and resuming an estimation run.

## Run the checks

From the package root, open MATLAB and run:

```matlab
cd replication
validate_package('saved')
validate_package('smoke')
```

The saved check verifies the package snapshot, regenerates all figures and tables, checks their exports and compares selected numerical values with the paper. The smoke check stages a fresh estimation workspace, reconstructs the simulation designs, fits selected empirical models with three bootstrap draws, and runs one replication of each parametric reference design.

For an independent numerical audit, run the following from the package root with Python, NumPy and pandas installed:

```sh
python validation/check_numerics.py
```

The audit independently refits the quadratic, changepoint and local-quadratic estimators to the frozen monthly data. It also checks fitted paths, crossing dates and archived simulation summaries.

For file-hashing, corruption-detection and estimation-resume checks, run:

```sh
matlab -nojvm -batch "addpath('validation'); check_java_free"
```

This command tests execution without a Java runtime. It checks known SHA-256 digests and creates a separate validation run under `outputs/reestimated/`. It fits the ERA5 quadratic and local-quadratic models, resumes the same run and tests rejection of altered input and source files in disposable copies.

## Numerical validation results

| Check | Result |
|---|---|
| Figures and tables | All 36 entry points correspond to the main paper and supplement. The saved workflow exports all 13 figures and 23 tables and passes the numerical checks. For Tables S7 and S8, all 96 displayed coverage, Monte Carlo error and tail-probability cells, including structural dashes, agree with the paper to its printed precision. |
| Independent empirical refits | The quadratic trend, grid-searched single changepoint and partially linear local-quadratic estimator agree with the archived monthly estimates within an absolute tolerance of 1e−8 for all four reported quantities, all five records and their separately fitted mean. All six selected breakpoint dates agree. |
| Fitted paths and sensitivity | All 24 sets of endpoint levels, rates, accelerations and past-decade rate differences agree with their saved fit paths. The method ordering and the lower-order smoother's five within-record comparisons agree with the paper, including its mean rate of 0.389 °C per decade. |
| Crossing dates | All twelve conditional 1.5 °C and 2 °C crossing dates, including the flexible constant-rate alternatives, agree with Tables S2 and S3 to their printed precision when calculated from unrounded mean-record fits. |
| Simulation summaries | Two-sided coverage, Monte Carlo standard errors and lower/upper miss frequencies in Tables S7 and S8 agree with calculations from 9,600 estimator–estimand–replication rows across 1,200 simulated records. The structural changepoint dash, all 80 accuracy summary rows and the 16/20 RMSE ranking also agree. |
| State-space examples | Annual ERA5 point estimates agree with the archive within 1e−6, and all three test bootstrap draws succeed. A reduced monthly ERA5 example without covariates matches the reference smoothed state paths and three-draw bootstrap bands exactly under the same ARMA(1,1) specification, settings and seed. |
| Estimation and resume | The ERA5 quadratic and local-quadratic example reproduces all eight reported point estimates within 2.19e−15 of the archive. Resuming the same RunID reuses both completed fits with identical summaries and cache-file digests. |

## File integrity

The package snapshot lists 178 files in `data/SNAPSHOT.csv`. The estimation-workspace manifest lists 61 files in `replication/engines/SHA256SUMS.tsv`. Both manifests pass independent verification. File contents extracted from the distribution ZIP match the snapshot.

The native MATLAB SHA-256 implementation passes 20 independent expected-digest checks: empty input, standard text vectors, one million `a` bytes, and binary inputs spanning padding, block and 1 MiB read-chunk boundaries. Expected digests in `validation/sha256_vectors.csv` were computed independently with Python's `hashlib.sha256`. Filenames containing spaces, quotes and non-ASCII characters also pass.

The snapshot verifier rejects a same-size one-byte change. Estimation-workspace verification rejects same-size changes to both an input file and a MATLAB source file. Missing-file hashing fails as expected. Resuming a run requires the same package, inputs and RunID; different inputs or code require a new RunID.

The bundled hashing core's source attribution and licence are in [DataHash_NOTICE.md](replication/engines/methods/private/DataHash_NOTICE.md) and [DataHash_LICENSE.txt](replication/engines/methods/private/DataHash_LICENSE.txt).

## Tested environment and scope

The saved and smoke checks pass in MATLAB R2025b. The saved workflow and the file-integrity/estimation-resume check also pass on macOS (Apple silicon) with MATLAB R2025b Update 4 started with `-nojvm`, with `usejava('jvm') == 0`. Other MATLAB releases and operating systems are not covered by these execution checks.

Fresh-estimation checks use small bootstrap samples. They do not constitute a complete rerun of the 1,999-draw empirical fits or the production Monte Carlo experiments. The supplied state-space fits, uncertainty intervals and Monte Carlo results are archived production output. The independent refits cover the other three estimators' monthly point estimates.
