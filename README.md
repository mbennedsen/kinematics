# Conventional trend models may understate the current rate of global warming

MATLAB replication package by Mikkel Bennedsen and J. Eduardo Vera-Valdés. The monthly vintage ends in August 2026, with July endpoints for HadCRUT5 and the separately fitted arithmetic-mean series.

## Reproduce the paper from saved estimates

Extract the ZIP or clone the repository, then open MATLAB in the folder containing this README. Run

```matlab
cd replication
verify_snapshot
figure3('saved')
TableS1('saved')
run_all('saved')
```

Each numbered entry point defaults to saved mode. The scripts build figures and tables from numerical results. Find PDF, PNG and editable FIG figures and HTML, CSV, MAT and LaTeX tables under `outputs/reproduced/saved`. MATLAB supplies the typography. To keep a figure open, use, for example, `figure3('saved',struct('Visible','on'))`. The saved workflow uses base MATLAB with graphics export support. The package is tested with MATLAB R2025b Update 4.

## Re-estimate from the frozen data

```matlab
% From replication/
runOptions = struct('Workers',1,'RunID','paper_run');
TableS4('reestimate',runOptions)  % Regenerate simulation truths first
figure3('reestimate',runOptions)  % Fit the four models to the arithmetic mean
run_all('reestimate',runOptions)  % Production workload
```

Fresh mode reads `data/` and the bundled engines. It writes its own fits, checkpoints, tables and figures to `outputs/reestimated/<RunID>/`. Reuse an ID for successive calls in the same run. Give concurrent MATLAB sessions distinct IDs. Omit the ID for a new run. Allow write access and sufficient disk space for bootstrap checkpoints. Re-estimation requires Optimization Toolbox and Statistics and Machine Learning Toolbox. Multiple workers also require Parallel Computing Toolbox.

Empirical bootstrap intervals use 1,999 draws where specified. Accuracy uses five designs with 1,000 replications each. S7 uses 400 outer records of the same smooth path for all four methods, with 399 inner draws for bootstrap intervals. S8 uses 400 records per matching parametric path and 399 inner draws for the changepoint intervals. The quadratic intervals use HAC inference. These are substantial computations. A fresh fit can differ slightly with the optimizer, MATLAB release and worker count.

For S7 and S8 the package includes the replication-level results as well as summaries. S7 uses the fixed integrated NASA path from manuscript DGP 3. S8 uses the quadratic and February 2014 broken-linear paths and the matching benchmark for each. Unused third-party columns are marked missing, as detailed in [DATA_SOURCES.md](DATA_SOURCES.md).

## Package layout

| Location | Contents |
|---|---|
| `replication/` | Numbered MATLAB entry points, `run_all` and the validation commands |
| `replication/engines/methods/` | Shared estimators, bootstrap routines and simulation path generators |
| `replication/engines/empirical/` and `sensitivity/` | Monthly and annual fits, the no-covariate exercise and sensitivity analyses |
| `replication/engines/accuracy/`, `smooth_coverage/` and `reference_coverage/` | Experiments for Tables S5–S6, S7 and S8, respectively |
| `data/` | Frozen analysis panels and model settings, including the shared simulation calibration in `Simulation/` |
| `outputs/saved/` | Numerical results used to render the published figures and tables |

Each re-estimation run verifies its input and source checksums before fitting. To resume, use the same package and RunID. Use a new RunID if the inputs or code differ.

## Manuscript output map

| Entry point | Manuscript output |
|---|---|
| figure1 | Warming level, rate and acceleration for the mean record since 2000 |
| figure2 | Endpoint estimates and intervals by method and record |
| figure3 | Conditional continuations to 2 °C with a 1.5 °C inset |
| FigureS1, FigureS2 | Temperature observations and climate covariates |
| FigureS3–FigureS8 | Full-sample paths for mean, NASA, NOAA, HadCRUT5, Berkeley Earth and ERA5, in order |
| FigureS9 | Continuations including the two flexible constant-rate alternatives |
| FigureS10 | Simulation acceleration paths |
| TableS1 | Monthly endpoint and past-decade estimates |
| TableS2, TableS3 | Conditional 1.5 °C and 2 °C crossing dates |
| TableS4–TableS6 | Simulation truth, RMSE and bias |
| TableS7, TableS8 | Coverage under the smooth NASA path and matching parametric paths |
| TableS9 | Annual-data estimates |
| TableS10, TableS11 | Quadratic and changepoint starting-year sensitivity |
| TableS12, TableS13 | Nonparametric bandwidth and state-space smoothing sensitivity |
| TableS14 | Lower-order spline smoother |
| TableS15 | Alternative starting dates for average acceleration |
| TableS16 | Semiparametric estimates without climate covariates |
| TableS17, TableS18 | Dynamic parameters and covariate basis coefficients |
| TableS19, TableS20 | Monthly variance scales and residual innovation variances |
| TableS21–TableS23 | Lag weights, standardized-innovation diagnostics and changepoint-bootstrap diagnostics |

## Data, citation and verification

The primary inputs are the frozen monthly and annual analysis panels. The mean record is constructed after rebasing each product to a 1991–2020 mean of 0.88 °C. Monthly past-decade acceleration uses exactly 120 lags. The annual exercise ends in 2025. The simulation calibrations and annual smoothing parameters are fixed as in the manuscript. This package starts at the assembled analysis data. Provider downloads and upstream construction are described in [DATA_SOURCES.md](DATA_SOURCES.md).

See [VALIDATION.md](VALIDATION.md) for the validation results and testing scope. From the replication folder, `validate_package('saved')` checks file integrity, regenerates all 36 items and verifies exported files and numerical anchors. `validate_package('smoke')` exercises fresh estimation and simulation design. The numbered workflows run locally in MATLAB using the bundled inputs. An optional independent audit using Python, NumPy and pandas is in `validation/check_numerics.py`.

## Licence and citation

The authors’ MATLAB and Python code and accompanying software documentation are available under the [MIT License](LICENSE). Keep the copyright and permission notices with copies or substantial portions of the software. The bundled DataHash SHA-256 core has separate [BSD 3-Clause licence terms](replication/engines/methods/private/DataHash_LICENSE.txt), with source attribution in [DataHash_NOTICE.md](replication/engines/methods/private/DataHash_NOTICE.md). Third-party observations in the panels and MATLAB calibrations retain the [data providers’ terms](DATA_SOURCES.md).

If you use the methods, code or results in research, please cite the paper. The public preprint is cited under its arXiv title, *The kinematics of global warming: Semiparametric analysis of warming level, rate, and acceleration*. [CITATION.cff](CITATION.cff) supplies this preferred citation.

Bennedsen, Mikkel, and J. Eduardo Vera-Valdés. 2026. “The kinematics of global warming: Semiparametric analysis of warming level, rate, and acceleration.” arXiv:2608.18298. [https://doi.org/10.48550/arXiv.2608.18298](https://doi.org/10.48550/arXiv.2608.18298)
