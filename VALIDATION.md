# Validation of the v60 replication package

Assessment on 20 November 2026. The public interface comprises 3 main figures, 10 supplementary figures and 23 supplementary tables.

| Check | Evidence |
|---|---|
| Numbering and source data | All 36 entry points mapped to labels in the two compiled v60 manuscripts. The numeric table bodies of 21 supplementary tables match v57. The two changed tables, S7 and S8, were checked against their respective v60 experiments. All 96 displayed coverage, Monte Carlo error and tail-probability cells, including structural dashes, agree to manuscript rounding. |
| Independent empirical refits | The quadratic trend, the grid-searched single changepoint and the partially linear local-quadratic estimator were fitted independently in Python to the frozen monthly inputs. All four reported quantities matched the archived estimates for each of the five records and their separately fitted mean within an absolute tolerance of 1e−8. The six selected breakpoint dates also matched. |
| Path and sensitivity checks | Checked all 24 sets of endpoint levels, rates, accelerations and past-decade rate differences against their saved fit paths. Confirmed the method ordering and the lower-order smoother’s five within-record comparisons, as well as its mean rate of 0.389 °C per decade. |
| Crossing dates | Recomputed the twelve 1.5 °C and 2 °C crossing dates, including both flexible constant-rate alternatives, from the unrounded mean-record fits. Results agree with Tables S2/S3 to their printed precision. |
| Simulations | Reconstructed two-sided coverage, its Monte Carlo standard error and lower/upper miss frequencies for S7 and S8 from 9,600 estimator–estimand–replication rows across 1,200 simulated records. Verified the structural changepoint dash. Checked all 80 archived accuracy summary rows and the 16/20 RMSE ranking. |
| Integrity and portability | SHA-256 snapshot for distributed files. Clean ZIP extraction, archived-file verification and independent numerical audit. The published engines’ input and source checks were independently inspected in their staged layout. The portable calibration uses a relative output directory. Six unused observational vectors are explicitly missing in all three distributed copies of the MATLAB calibration. The MATLAB header, object structure, references, attributes and every other dataset were checked against the previous calibration. The retained CSV fields were checked against their original text at every row. |
| MATLAB code | The MATLAB R2025b Code Analyzer parsed the sources. Mikkel Bennedsen also ran the two v60 validation commands in a local MATLAB session and reported both passes on 20 November 2026. |

## Native MATLAB validation

The author reported `VALIDATION PASS: saved` and `VALIDATION PASS: smoke`. The saved check regenerated 13 figures and 23 tables, verified the exported files, and checked selected crossing dates and coverage estimates. The smoke check exercised selected fresh empirical fits with reduced bootstrap draws, checked the simulation designs, and ran one outer replication for each of the two parametric reference designs. A full production rerun and independent coauthor validation are pending. The local MATLAB release was not reported.

After this run, Figure 1b was updated to draw the changepoint rate as a step at the fitted breakpoint. This plotting change is included in the next independent MATLAB check.

To repeat these checks from the extracted package:

```matlab
cd replication
validate_package('saved')
validate_package('smoke')
```

For the independent audit, run `python validation/check_numerics.py` from the package root with NumPy and pandas installed. The saved state-space fits, intervals and Monte Carlo results are archived production output. The independent refits cover the other three estimators’ monthly point estimates.

## Public-deposit assessment

The package is ready for independent review from a fresh repository clone. The authors’ code is released under the [MIT License](LICENSE). The [data acknowledgements](DATA_SOURCES.md) identify the products, transformations and applicable third-party terms. When the public repository issues a persistent identifier, add it to the software citation metadata and the manuscript’s availability statements.
