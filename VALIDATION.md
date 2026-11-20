# Validation of the v60 replication package

Assessment on 20 November 2026. The public interface comprises 3 main figures, 10 supplementary figures and 23 supplementary tables.

| Check | Evidence |
|---|---|
| Numbering and source data | All 36 entry points mapped to labels in the two compiled v60 manuscripts. For Tables S7 and S8, the replication records were checked against their respective experiments. All 96 displayed coverage, Monte Carlo error and tail-probability cells, including structural dashes, agree to manuscript rounding. |
| Independent empirical refits | The quadratic trend, the grid-searched single changepoint and the partially linear local-quadratic estimator were fitted independently in Python to the frozen monthly inputs. All four reported quantities matched the archived estimates for each of the five records and their separately fitted mean within an absolute tolerance of 1e−8. The six selected breakpoint dates also matched. |
| Path and sensitivity checks | Checked all 24 sets of endpoint levels, rates, accelerations and past-decade rate differences against their saved fit paths. Confirmed the method ordering and the lower-order smoother’s five within-record comparisons, as well as its mean rate of 0.389 °C per decade. |
| Crossing dates | Recomputed the twelve 1.5 °C and 2 °C crossing dates, including both flexible constant-rate alternatives, from the unrounded mean-record fits. Results agree with Tables S2/S3 to their printed precision. |
| Simulations | Reconstructed two-sided coverage, its Monte Carlo standard error and lower/upper miss frequencies for S7 and S8 from 9,600 estimator–estimand–replication rows across 1,200 simulated records. Verified the structural changepoint dash. Checked all 80 archived accuracy summary rows and the 16/20 RMSE ranking. |
| Integrity and portability | The release ZIP was extracted and checked from a fresh local Git clone. The file manifest and the estimation-workspace manifest were verified independently. Every frozen analysis input and calibration matches the audited snapshot byte for byte. Two saved CSVs omit internal source-file metadata. Every retained field matches its original text in every row. |
| MATLAB code | MATLAB R2025b Code Analyzer parsed the distributed sources. |

## Native MATLAB validation

The saved check regenerates all 13 figures and 23 tables, checks their exports and verifies numerical anchors. The smoke check fits selected empirical models with reduced bootstrap draws, reconstructs the simulation designs, and runs one replication for each parametric reference design. Both commands passed in MATLAB R2025b on 20 November 2026 from a clean local Git clone of the reorganised package.

Two additional checks exercised the state-space code. For annual ERA5 data, the refitted point estimates agreed with the archived values to within 10⁻⁶ and all three test bootstrap draws succeeded. For the monthly fit without covariates, a reduced ERA5 example used the same ARMA(1,1) specification, settings and random seed in both implementations. The smoothed state paths and the three-draw bootstrap path bands matched exactly.

These checks use small bootstrap samples. The full 1,999-draw empirical fits and production Monte Carlo experiments have not been rerun for this release.

To repeat these checks from the extracted package:

```matlab
cd replication
validate_package('saved')
validate_package('smoke')
```

For the independent audit, run `python validation/check_numerics.py` from the package root with NumPy and pandas installed. The saved state-space fits, intervals and Monte Carlo results are archived production output. The independent refits cover the other three estimators’ monthly point estimates.

## Public-deposit assessment

The package is ready for independent review from a fresh repository clone. The authors’ code is released under the [MIT License](LICENSE). The [data acknowledgements](DATA_SOURCES.md) identify the products, transformations and applicable third-party terms. When the public repository issues a persistent identifier, add it to the software citation metadata and the manuscript’s availability statements.
