
# ctrest

**ctrest** is an R package for estimating the density of ground-dwelling
mammals and birds from **camera-trap video data** using the **REST**
(Random Encounter and Staying Time) model or the **RAD-REST** (Random
Encounter and Staying Time model Relying on All Detections) model.

Both models are non-invasive, require no individual recognition, and
work with standard video recordings from fixed camera traps.

------------------------------------------------------------------------

## Background

The REST model estimates animal density from two key measurements
recorded for every video clip:

1.  **Number of passes** — how many times the animal entered the focal
    area (the area clearly visible in the camera’s field of view) in
    that clip.
2.  **Staying time** — how long the animal spent inside the focal area.

Combined with **camera trapping effort** (days) and an estimate of the
**proportion of time the animal is active**, the REST model derives
density without requiring a conventional distance-sampling survey
(Nakashima et al. 2018).

The **RAD-REST** model is an extension that relaxes the requirement of
counting passes in every video: it is sufficient to count passes in a
random subsample, which greatly reduces the annotation workload
(Nakashima et al. 2026).

------------------------------------------------------------------------

## Installation

``` r
# install.packages("devtools")
devtools::install_github("YoshihiroNakashima/ctrest", build_vignettes = TRUE)
```

For a guided tour after installation:

``` r
library(ctrest)
vignette("ctrest")
#> Warning: vignette 'ctrest' not found
```

### System requirements

- R ≥ 4.1
- [NIMBLE](https://r-nimble.org/) is used internally for Bayesian MCMC.
  On Windows, [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  must be installed. On macOS, Xcode Command Line Tools are required.
- A multi-core CPU is strongly recommended because chains are run in
  parallel.

------------------------------------------------------------------------

## Choosing between REST and RAD-REST

| Situation                                      | Recommended model |
|------------------------------------------------|-------------------|
| Passes are counted for every video             | REST              |
| Passes are counted for a random subsample only | RAD-REST          |
| You want to minimise annotation time           | RAD-REST          |
| You want the simplest possible workflow        | REST              |

REST data can always be analysed with RAD-REST (but not vice versa).

------------------------------------------------------------------------

## Workflow overview

    detection_data ──┬── format_station_data() ──── add_effort() ──────────────────┐
                     ├── format_stay()                                               │
                     │         └── bayes_stay_selection()   (optional)              ├──▶ bayes_rest()
                     └── format_activity()                                           │
    station_data   ──┘                                                               │
                                                                                     └──▶ bayes_rest_multi()

The analysis has five stages:

1.  **Format station data** (`format_station_data` + `add_effort`)
2.  **Format staying time data** (`format_stay`)
3.  **Format activity data** (`format_activity`)
4.  *(Optional)* **Select the staying time distribution**
    (`bayes_stay_selection`)
5.  **Estimate density** (`bayes_rest` or `bayes_rest_multi`)

------------------------------------------------------------------------

## Example datasets

Two example datasets ship with the package:

``` r
library(ctrest)

head(detection_data)
head(station_data)
```

**`detection_data`** has one row per video clip and contains:

| Column | Description |
|----|----|
| `Station` | Camera station ID |
| `DateTime` | Date and time of detection |
| `Term` | Survey round / visit ID |
| `Species` | Species detected |
| `y` | Number of passes through the focal area in this clip |
| `Stay` | Staying time (seconds) inside the focal area |
| `Cens` | Censoring flag: `1` = right-censored (animal still in frame at clip end), `0` = observed |

**`station_data`** has one row per camera station and contains:

| Column    | Description                                          |
|-----------|------------------------------------------------------|
| `Station` | Camera station ID (matches `detection_data$Station`) |
| `x1`      | A continuous environmental covariate                 |
| `x2`      | A categorical environmental covariate                |

------------------------------------------------------------------------

## Stage 1: Format station data

`format_station_data()` aggregates the per-video pass counts into
**station-level summaries** and joins station covariates.

For the REST model it produces a column `Y` (total passes per station).
For RAD-REST it additionally produces columns `y_0`, `y_1`, `y_2`, …
(counts of videos with 0, 1, 2, … passes observed).

``` r
# REST
station_data_rest <- format_station_data(
  detection_data   = detection_data,
  station_data     = station_data,
  col_name_station = "Station",
  col_name_species = "Species",
  col_name_y       = "y",
  model            = "REST"
)

head(station_data_rest)
```

``` r
# RAD-REST
station_data_rad <- format_station_data(
  detection_data   = detection_data,
  station_data     = station_data,
  col_name_station = "Station",
  col_name_species = "Species",
  col_name_y       = "y",
  model            = "RAD-REST"
)
#> 400 record(s) with NA in 'y' were removed before aggregation.

head(station_data_rad)
#> # A tibble: 6 × 9
#>   Station Species     N   y_0   y_1   y_2   y_3     x1 x2   
#>   <chr>   <chr>   <int> <int> <int> <int> <int>  <dbl> <chr>
#> 1 ST001   SP01        5     1     4     0     0  0.520 C    
#> 2 ST002   SP01        1     1     0     0     0  0.445 C    
#> 3 ST003   SP01       11     6     4     1     0 -1.17  A    
#> 4 ST004   SP01       17     5    12     0     0  0.319 A    
#> 5 ST005   SP01        4     1     1     1     1  0.332 C    
#> 6 ST006   SP01        7     2     3     1     1 -2.02  B
```

### Compute camera trapping effort

`add_effort()` calculates **how many days each camera was operating** by
computing the span between the first and last detection at each station
(within each survey term, if `col_name_term` is supplied).

Set `plot = TRUE` to display a Gantt-style timeline of camera operation
— a useful sanity check before running models.

``` r
station_effort_rest <- add_effort(
  detection_data         = detection_data,
  station_data_formatted = station_data_rest,
  col_name_station       = "Station",
  col_name_datetime      = "DateTime",
  col_name_term          = "Term",   # NULL if there are no survey rounds
  plot                   = TRUE
)

head(station_effort_rest)
```

> **Note.** The effort estimate (last detection − first detection) is an
> approximation. For more accurate estimates, provide actual camera
> operation logs directly.

------------------------------------------------------------------------

## Stage 2: Format staying time data

`format_stay()` joins the detection records to station metadata, renames
the key columns to the standard names expected by the model functions,
and flags any data quality issues.

``` r
stay_data <- format_stay(
  detection_data   = detection_data,
  station_data     = station_data,
  col_name_station = "Station",
  col_name_species = "Species",
  col_name_stay    = "Stay",
  col_name_cens    = "Cens"
)
#> 400 record(s) with missing values in 'Stay' or 'Cens' were excluded.

head(stay_data)
#> # A tibble: 6 × 6
#>   Station Species  Stay  Cens    x1 x2   
#>   <chr>   <chr>   <dbl> <int> <dbl> <chr>
#> 1 ST001   SP01      6.7     0 0.520 C    
#> 2 ST001   SP01      9.6     0 0.520 C    
#> 3 ST001   SP01      2.1     0 0.520 C    
#> 4 ST001   SP01      4.4     0 0.520 C    
#> 5 ST001   SP01      3.3     0 0.520 C    
#> 6 ST002   SP01      2.4     0 0.445 C
```

The output contains `Station`, `Species`, `Stay` (seconds), `Cens`
(0/1), and all covariates from `station_data`. Rows with missing or
non-positive staying times are removed automatically, with a message
reporting how many were dropped.

------------------------------------------------------------------------

## Stage 3: Format activity data

`format_activity()` converts detection datetimes to **radians on a \[0,
2π\] scale** so that 00:00 = 0 and 23:59 ≈ 2π. It also filters for
independent detections (one detection per station per species per
`indep_time`-minute window).

``` r
activity_data <- format_activity(
  detection_data    = detection_data,
  col_name_station  = "Station",
  col_name_species  = "Species",
  col_name_datetime = "DateTime",
  indep_time        = 30   # minimum minutes between independent detections
)

head(activity_data)
#> # A tibble: 6 × 3
#>   Species  Station  time
#>   <chr>    <chr>   <dbl>
#> 1 Surveyor ST001    2.55
#> 2 SP11     ST001    2.90
#> 3 SP03     ST001    2.90
#> 4 SP11     ST001    3.22
#> 5 SP11     ST001    1.03
#> 6 SP03     ST001    3.22
```

------------------------------------------------------------------------

## Stage 4 (optional): Select the staying time covariate model

`bayes_stay_selection()` compares candidate **covariate models** for
staying time via WAIC for a **single** chosen probability distribution.
To also compare distributions, call the function separately for each
candidate (e.g. `"lognormal"`, `"gamma"`, `"weibull"`) and compare the
best WAIC values across runs.

Supported distributions for `stay_family`: `"exponential"`, `"gamma"`,
`"weibull"`, `"lognormal"`.

``` r
stay_sel <- bayes_stay_selection(
  formula_stay       = Stay ~ 1 + x1,   # covariates on mean staying time
  random_effect_stay = NULL,             # or "Station" for a station-level random intercept
  stay_data          = stay_data,
  col_name_cens      = "Cens",
  stay_family        = "lognormal",      # distribution to evaluate
  iter               = 5000,
  warmup             = 1000,
  chains             = 3,
  thin               = 4,
  all_comb           = TRUE,             # if TRUE, compare all covariate subsets of formula_stay
  target_species     = "SP01"
)

print(stay_sel)
```

The printed output shows:

- A WAIC table (columns: `Model`, `Family`, `Random_effect`, `WAIC`,
  `Note`) comparing all covariate subsets for the chosen distribution,
  with the best model marked.
- Posterior estimates of the mean staying time from the best model.
- Convergence diagnostics (Rhat).
- A Bayesian p-value assessing model fit.

You can access the raw MCMC samples for further diagnostics:

``` r
MCMCvis::MCMCtrace(stay_sel$samples, pdf = FALSE)
```

------------------------------------------------------------------------

## Stage 5: Estimate density — REST model

Pass the data prepared in Stages 1–4 to `bayes_rest()` with
`model = "REST"`. Set `focal_area` to the size of the camera’s focal
area in **square metres**.

``` r
result_rest <- bayes_rest(
  formula_stay        = Stay ~ 1 + x1,
  formula_density     = ~ 1,             # density model (no LHS)
  formula_enter       = ~ 1,             # required argument; ignored for REST
  station_effort_data = station_effort_rest,
  stay_data           = stay_data,
  activity_data       = activity_data,
  random_effect_stay  = NULL,
  activity_estimation = "kernel",        # "kernel" or "mixture"
  bw_adj              = 1.0,             # bandwidth adjustment for kernel estimation
  stay_family         = "lognormal",     # from bayes_stay_selection output
  focal_area          = 1.96,            # focal area in m²
  iter                = 5000,
  warmup              = 1000,
  chains              = 3,
  thin                = 2,
  model               = "REST",
  all_comb            = FALSE,           # TRUE = compare all density covariate subsets
  target_species      = "SP01"
)
```

### Activity estimation: kernel vs mixture

The `activity_estimation` argument controls how the proportion of time
the animal is active is estimated:

| Value | Method | Notes |
|----|----|----|
| `"kernel"` | Fixed kernel density (Rowcliffe et al. 2014) | Fast; `bw_adj` controls bandwidth |
| `"mixture"` | Nonparametric von Mises mixture (Nakashima et al. 2025) | Fully Bayesian; propagates uncertainty; `C` sets max components |

The mixture method is recommended when activity patterns are complex
(e.g., multimodal) or when you want activity uncertainty to propagate
into the density estimate.

### Interpreting the output

``` r
print(result_rest)
```

The printed summary shows:

    === ctrest: Density estimation ===
    Species    : SP01
    Model      : REST
    Stay family: lognormal

    --- Model comparison (WAIC) ---
      Model  random_effect_stay   WAIC  Note
      ~ 1    NULL               1234.5  <- best

    --- Posterior estimates (best model) ---
      Species  Station  Variable    mean    sd   lower  median   upper  Rhat  n.eff     cv
      SP01     All      density    0.052  0.008  0.038   0.051   0.069  1.00   3000  0.154
      SP01     All      mean_stay  4.230  0.312  3.640   4.218   4.872  1.00   3000  0.074

    --- Convergence ---
      All 12 monitored parameter(s): Rhat <= 1.1.

    Note: Full MCMC samples : $samples
          Long-format samples: $tidy_samples
          ...

- **`density`** (individuals per km²) — divide by 100 to convert to
  individuals per ha.
- **`mean_stay`** — posterior mean of mean staying time in seconds.
- **Rhat** — the Gelman-Rubin convergence diagnostic. Values ≤ 1.1
  indicate good convergence. If Rhat is large, increase `iter` or
  `chains`.
- **n.eff** — effective sample size. Values ≥ 1000 are generally
  acceptable.
- **cv** — coefficient of variation (sd / mean).

### Accessing raw MCMC samples

``` r
# Full posterior samples (coda::mcmc.list)
result_rest$samples

# Long-format tibble — useful for ggplot2
result_rest$tidy_samples

# Trace plots
MCMCvis::MCMCtrace(result_rest$samples, pdf = FALSE)

# Density plot of a parameter
MCMCvis::MCMCplot(result_rest$samples, params = "density")

# Full posterior summary table
print(result_rest$summary_result)
```

### Activity curve (mixture method only)

When `activity_estimation = "mixture"`, the output includes
`$activity_curve`, a data frame of the posterior activity density
function evaluated on a fine time grid. This can be plotted directly
(replace `result_mix` with your own mixture-based result):

``` r
library(ggplot2)

ggplot(result_mix$activity_curve, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = mean), colour = "steelblue") +
  scale_x_continuous(
    name   = "Time of day",
    breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
    labels = c("00:00", "06:00", "12:00", "18:00", "24:00")
  ) +
  ylab("Activity density") +
  theme_bw()
```

### Covariates in the density model

When covariates are included in `formula_density` (or `formula_stay`),
one row per station is returned in `$summary_result`, showing
station-specific estimates:

``` r
result_cov <- bayes_rest(
  formula_stay        = Stay ~ 1,
  formula_density     = ~ 1 + x1,        # density varies with covariate x1
  formula_enter       = ~ 1,
  station_effort_data = station_effort_rest,
  stay_data           = stay_data,
  activity_data       = activity_data,
  activity_estimation = "kernel",
  stay_family         = "lognormal",
  focal_area          = 1.96,
  iter                = 5000,
  warmup              = 1000,
  chains              = 3,
  thin                = 2,
  model               = "REST",
  target_species      = "SP01"
)

print(result_cov)
# One row per station
print(result_cov$summary_result)
```

Use `all_comb = TRUE` to automatically compare all possible covariate
subsets via WAIC and select the best model:

``` r
result_allcomb <- bayes_rest(
  formula_stay        = Stay ~ 1,
  formula_density     = ~ 1 + x1 + x2,   # candidate covariates
  formula_enter       = ~ 1,
  station_effort_data = station_effort_rest,
  stay_data           = stay_data,
  activity_data       = activity_data,
  activity_estimation = "kernel",
  stay_family         = "lognormal",
  focal_area          = 1.96,
  iter                = 5000,
  warmup              = 1000,
  chains              = 3,
  thin                = 2,
  model               = "REST",
  all_comb            = TRUE,             # fit all 4 subsets: ~1, ~x1, ~x2, ~x1+x2
  target_species      = "SP01"
)
print(result_allcomb)
# WAIC table for all candidate models
print(result_allcomb$WAIC)
```

------------------------------------------------------------------------

## Stage 5: Estimate density — RAD-REST model

The RAD-REST model extends REST by explicitly modelling the number of
animal passes **per video clip** using a count distribution. This
requires station data formatted with `model = "RAD-REST"` (which retains
per-video counts) and adds `formula_enter` to specify covariates on the
mean number of passes per video.

First, compute camera trapping effort from the RAD-REST formatted
station data:

``` r
station_effort_rad <- add_effort(
  detection_data         = detection_data,
  station_data_formatted = station_data_rad,
  col_name_station       = "Station",
  col_name_datetime      = "DateTime",
  col_name_term          = "Term"
)
```

Then run the model:

``` r
result_rad <- bayes_rest(
  formula_stay        = Stay ~ 1,
  formula_density     = ~ 1,
  formula_enter       = ~ 1,             # covariates on mean passes per video
  station_effort_data = station_effort_rad,
  stay_data           = stay_data,
  activity_data       = activity_data,
  activity_estimation = "kernel",
  stay_family         = "lognormal",
  focal_area          = 1.96,
  iter                = 5000,
  warmup              = 1000,
  chains              = 3,
  thin                = 2,
  model               = "RAD-REST",
  all_comb            = FALSE,
  target_species      = "SP01"
)
```

### Interpreting the output

``` r
print(result_rad)
```

The output structure is identical to the REST model. The key addition
is:

- **`mean_pass`** — posterior mean of the mean number of passes per
  video clip (RAD-REST only).
- **`density`** (individuals per km²) — computed via the same formula as
  REST; units and interpretation are unchanged.

------------------------------------------------------------------------

## Standalone activity estimation

You can also estimate activity patterns independently of density using
`bayes_activity()`. This is useful for comparing activity patterns
across species or seasons, or for quick checks before running the full
model.

``` r
act_result <- bayes_activity(
  activity_data  = activity_data,
  C              = 10,       # maximum number of von Mises components
  iter           = 5000,
  warmup         = 1000,
  chains         = 3,
  thin           = 2,
  target_species = "SP01"
)
#> Compiling the model. This may take a moment...
#> Running MCMC sampling. Please wait...
#> Estimation is finished!
```

Print a text summary (no plot):

``` r
printResultActivity(act_result, plot = FALSE, bw_adj = 1.0)
#> 
#> === ctrest: Activity estimation ===
#> Species: SP01
#> 
#> --- Activity proportion ---
#> # A tibble: 1 × 9
#>   Variable             mean     sd lower median upper  Rhat n.eff     cv
#>   <chr>               <dbl>  <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>  <dbl>
#> 1 activity_proportion 0.383 0.0151 0.355  0.383 0.415  1.01  3142 0.0394
#> 
#> --- Convergence ---
#>   All 1 monitored parameter(s): Rhat <= 1.1.
#> 
#> --- Bayesian p-value ---
#>   1.000  (model fit appears poor)
#> 
#> Note: Full MCMC samples : $samples
#>       Long-format samples: $tidy_samples
#>       Activity curve data: $activity_curve
#>       Trace plots        : MCMCvis::MCMCtrace(x$samples)
```

Print a summary with the activity curve plot:

``` r
printResultActivity(act_result, plot = TRUE, bw_adj = 1.0)
#> 
#> === ctrest: Activity estimation ===
#> Species: SP01
#> 
#> --- Activity proportion ---
#> # A tibble: 1 × 9
#>   Variable             mean     sd lower median upper  Rhat n.eff     cv
#>   <chr>               <dbl>  <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>  <dbl>
#> 1 activity_proportion 0.383 0.0151 0.355  0.383 0.415  1.01  3142 0.0394
#> 
#> --- Convergence ---
#>   All 1 monitored parameter(s): Rhat <= 1.1.
#> 
#> --- Bayesian p-value ---
#>   1.000  (model fit appears poor)
#> 
#> Note: Full MCMC samples : $samples
#>       Long-format samples: $tidy_samples
#>       Activity curve data: $activity_curve
#>       Trace plots        : MCMCvis::MCMCtrace(x$samples)
```

![](README_files/figure-gfm/print-activity-plot-1.png)<!-- -->

The plot shows: - **Blue ribbon**: 95% credible interval of the
posterior activity density. - **Blue line**: Posterior mean activity
density. - **Dashed grey line**: Non-parametric kernel density estimate
for comparison. - **Rug**: Raw detection times.

The activity proportion (probability that the animal is active at a
random moment) is reported in `$summary_result`.

------------------------------------------------------------------------

## Multi-species density estimation

`bayes_rest_multi()` estimates density for several species
simultaneously using a **hierarchical model** in which species-level
random effects are shared across species, stabilising estimates for rare
species.

``` r
result_multi <- bayes_rest_multi(
  formula_stay        = Stay ~ 1,
  formula_density     = ~ 1,
  formula_enter       = ~ 1,
  station_effort_data = station_effort_rad,   # RAD-REST format
  stay_data           = stay_data,
  activity_data       = activity_data,
  activity_estimation = "kernel",
  stay_family         = "lognormal",
  focal_area          = 1.96,
  iter                = 5000,
  warmup              = 1000,
  chains              = 3,
  thin                = 4,
  model               = "RAD-REST",
  target_species      = c("SP01", "SP02", "SP03")
)
#> Compiling the model. This may take a moment...
#> Running MCMC sampling. Please wait...
#> Estimation is finished!
print(result_multi)
#> 
#> === ctrest: Density estimation ===
#> Species    : SP01, SP02, SP03
#> Model      : RAD-REST
#> Stay family: lognormal
#> 
#> --- Model comparison (WAIC) ---
#>   WAIC: 12313.64
#> 
#> --- Posterior estimates ---
#> # A tibble: 9 × 11
#>   Species Station Variable    mean     sd lower median  upper  Rhat n.eff     cv
#>   <chr>   <chr>   <chr>      <dbl>  <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl>  <dbl>
#> 1 SP01    All     density[… 10.4   0.937  8.76  10.3   12.3    1.01  1849 0.0901
#> 2 SP02    All     density[…  6.94  0.648  5.78   6.91   8.37   1.01  1683 0.0934
#> 3 SP03    All     density[…  3.90  0.458  3.07   3.87   4.87   1     1449 0.117 
#> 4 SP01    All     mean_sta…  4.01  0.0875 3.84   4.01   4.18   1.02  2057 0.0218
#> 5 SP02    All     mean_sta…  4.20  0.121  3.97   4.19   4.45   1.02  1834 0.0288
#> 6 SP03    All     mean_sta…  3.97  0.138  3.71   3.96   4.25   1.01  1880 0.0347
#> 7 SP01    All     mean_pas…  0.921 0.0345 0.855  0.920  0.989  1     1959 0.0375
#> 8 SP02    All     mean_pas…  0.922 0.0386 0.845  0.921  0.999  1     2165 0.0419
#> 9 SP03    All     mean_pas…  0.942 0.0477 0.856  0.939  1.05   1     1040 0.0507
#> 
#> --- Convergence ---
#>   All 9 monitored parameter(s): Rhat <= 1.1.
#> 
#> Note: Full MCMC samples : $samples
#>       Long-format samples: $tidy_samples
#>       Trace plots        : MCMCvis::MCMCtrace(x$samples)
```

The output structure is identical to `bayes_rest()`. Density estimates
for all species are returned in `$summary_result`.

------------------------------------------------------------------------

## Video subsampling for RAD-REST

If you collected REST-format data (passes counted for all videos) but
want to apply RAD-REST to a subsample, use `select_videos()` to draw a
balanced random subsample:

``` r
sampled <- select_videos(
  detection_data    = detection_data,
  col_name_species  = "Species",
  col_name_station  = "Station",
  col_name_datetime = "DateTime",
  N_sampled         = 100,        # total videos to select
  Indep_criteria    = 30,         # independence interval in minutes
  target_species    = "SP01",
  seed              = 42
)
```

The function samples videos evenly across stations to avoid
over-representing highly active locations.

------------------------------------------------------------------------

## Function reference

| Function | Purpose |
|----|----|
| `format_station_data()` | Aggregate pass counts per station; join covariates |
| `add_effort()` | Compute camera trapping effort (days) per station |
| `format_stay()` | Prepare staying time data with censoring |
| `format_activity()` | Convert datetimes to radians; filter independent detections |
| `select_videos()` | Balanced random subsample for RAD-REST |
| `bayes_stay_selection()` | Compare staying time distributions via WAIC |
| `bayes_rest()` | Bayesian density estimation for a single species |
| `bayes_rest_multi()` | Hierarchical Bayesian density estimation for multiple species |
| `bayes_activity()` | Standalone Bayesian activity pattern estimation |
| `print.ResultStay()` | Print staying time model selection summary; called via `print()` |
| `print.ResultDensity()` | Print density estimation summary; called via `print()` |
| `print.ResultActivity()` | Print activity estimation summary without plot; called via `print()` |
| `printResultActivity()` | Print activity summary and optionally plot the activity curve (`plot = TRUE/FALSE`) |

------------------------------------------------------------------------

## Tips

**My Rhat values are \> 1.1.** Increase `iter` (e.g., to 10 000) and/or
`chains` (e.g., to 4). Also check whether the model is identifiable —
highly correlated parameters or very sparse data can hinder convergence.

**MCMC compilation takes a long time.** This is normal for NIMBLE; the
model is compiled once per chain. Subsequent runs with the same model
structure are not faster, but you can reduce `warmup` for exploratory
runs.

**The activity proportion estimate seems unreasonably high or low.**
Check `$activity_curve`: if the estimated density curve has unrealistic
peaks (e.g., very narrow spikes), try reducing `C` (fewer mixture
components) or using `activity_estimation = "kernel"` instead.

**Density estimates differ between REST and RAD-REST.** This is expected
when the subsample is small. RAD-REST inference is consistent but has
wider uncertainty intervals. With large subsamples the two should agree.

**I want to convert density to abundance.** Multiply `density`
(individuals / km²) by the area of your study region in km².

------------------------------------------------------------------------

## References

Nakashima, Y., Fukasawa, K., & Samejima, H. (2018). Estimating animal
density without individual recognition using information derivable
exclusively from camera traps. *Journal of Applied Ecology*, **55**(2),
900-910. <https://doi.org/10.1111/1365-2664.13059>

Nakashima, Y. Yajima, G. & Hongo, S. (2021). Estimating animal density
with camera traps: a practitioner’s guide of the REST model. *bioRxiv*.
<https://doi.org/10.1101/2021.05.18.444583>

Nakashima, Y., Yajima, G. & Matsuoka, R. (2026). Reducing data
processing effort in camera trap density estimation: Extending the REST
model by explicitly modelling animal detection processes. *Methods in
Ecology and Evolution*. <https://doi.org/10.1111/2041-210x.70248>

Rowcliffe, J. M., Kays, R., Kranstauber, B., Carbone, C., & Jansen, P.
A. (2014). Quantifying levels of animal activity using camera trap data.
*Methods in Ecology and Evolution*, **5**(11), 1170-1179.
<https://doi.org/10.1111/2041-210X.12278>
