# Drought, Wet Periods and Disinfection Byproduct Violations in Texas Drinking Water (2016-2025)

**Author:** Anamika Dristi, Ph.D.

**Status:** pilot analysis of public data, October 2026. It is not peer reviewed. The wet-period results were tested after the drought result and are exploratory.

## Question

When chlorine meets natural organic matter in source water, it forms disinfection byproducts (DBPs): total trihalomethanes (TTHM) and haloacetic acids (HAA5). These are regulated because long-term exposure is linked to bladder cancer and possibly to adverse reproductive outcomes. DBPs are the largest group of health-based drinking water violations in Texas and Louisiana (see the companion repository `drinking-water-violations-tx-la`).

My doctoral work showed that hydrologic extremes control how much organic carbon rivers carry. Drought cut carbon export from the Atchafalaya River, and treated wastewater adds dissolved organic carbon at two to five times river levels. This pilot asks whether those climate-driven changes in source water show up in public health records. Do DBP violations at Texas community water systems change with drought, or with wet periods that flush organic matter into rivers and reservoirs?

## Data (all public)

| Data | Source | Use |
|---|---|---|
| Texas community water systems and health-based DBP violations, 2015-2025 | EPA SDWIS via ECHO bulk download (processed in `drinking-water-violations-tx-la`) | 5,034 systems; TTHM and HAA5 MCL violations by start date |
| Weekly drought, percent of county area in each category | U.S. Drought Monitor county statistics (usdmdataservices.unl.edu) | County drought (D1 or worse), averaged over the past four quarters |
| Monthly county precipitation and temperature | NOAA nClimDiv | Quarterly precipitation anomaly (z-score against 1991-2020) and temperature anomaly |

## Methods

- **Panel:** each system in each quarter, 2016-2025 (about 179,000 system-quarters).
- **Outcome:** the number of TTHM or HAA5 MCL violations whose compliance period begins in that quarter.
- **Drought exposure:** county percent area in D1 or worse, averaged over the past four quarters. This matches the running-annual-average window that compliance is based on.
- **Wet exposures (exploratory):**
  - county precipitation anomaly in the current and previous quarter;
  - a "very wet" previous quarter (z > 1);
  - a "whiplash" term: a very wet quarter after a year with more than 30% of the county in drought.
- **Models:** Poisson regression with year, season, system size, source-water type and temperature anomaly as covariates, and standard errors clustered by county. Models are run for all systems, surface-water systems (own and purchased) and groundwater systems. A county fixed-effects version is included as a check.

## Results

**1. Surface water carries most of the burden.** Surface-water systems had 22.9 DBP MCL violations per 100 systems per year, compared with 6.3 for groundwater systems. That is about 3.7 times higher, consistent with more organic carbon in surface water.

**2. Drought itself was not associated with violations.** The rate ratio per 10 additional percentage points of county drought area was:

| Systems | Rate ratio | 95% CI |
|---|---|---|
| All | 0.97 | 0.88-1.07 |
| Surface water | 1.01 | 0.92-1.11 |
| Groundwater | 0.92 | 0.81-1.04 |
| All, with county fixed effects | 0.97 | 0.93-1.00 |

**3. Wet periods were associated with more violations (exploratory).**
- **All systems:** each one standard deviation of extra rainfall raised the violation rate by about 11%, in both the current quarter (1.12, 95% CI 1.04-1.19) and the previous quarter (1.11, 1.05-1.18). A very wet previous quarter raised it by about 31% (1.31, 1.12-1.53).
- **Surface-water systems:** a very wet quarter that *followed a drought year* added a further 61% (whiplash rate ratio 1.61, 1.11-2.33). This fits a flushing mechanism: organic matter builds up on dry land, then floods carry it into rivers and reservoirs.
- **Groundwater systems:** violations also rose after very wet quarters (1.60, 1.31-1.97). The reason is not clear from these data. Possible explanations include surface water reaching wells during floods, organic-rich aquifers, or patterns in monitoring and reporting.

**Figures:**
- `outputs/figures/drought_dbp_texas.png`: statewide quarterly DBP violations per 100 surface-water systems against drought area, with the drought rate ratios.
- `outputs/figures/dbp_drought_vs_wet_texas.png`: drought and wet-period rate ratios side by side.

## Limitations

- **Violations are not exposure.** They are compliance records, so they don't measure the DBP concentration people actually received, or any health outcome.
- **Exposures are measured at the county level.** County drought and rainfall stand in for conditions in each system's source watershed and reservoir. Systems that buy water may depend on rainfall far away.
- **Reporting can change over time.** Monitoring and enforcement practices vary by year, which the year terms only partly absorb.
- **Multiple testing.** The wet-period analyses were run after the drought result, so they are hypothesis-generating.

## Next steps

1. Link systems to their actual source watersheds and intakes, using USGS streamflow and TOC/DOC records where available.
2. Measure organic carbon flushing directly, using event-scale river DOC and satellite CDOM.
3. Add DBP concentration data (not only violations) and population served, to estimate exposure.
4. With mentorship in environmental epidemiology, link exposures to health records, for example county bladder cancer incidence or adverse birth outcomes. Ecological designs would be the first step.

## Run it

```bash
pip install pandas numpy statsmodels matplotlib
python scripts/01_analysis.py              # drought models, tables t1-t4, figure 1
python scripts/02_wet_periods_exploratory.py   # wet-period and whiplash models, table t5
python scripts/03_figures.py               # figure 2
```

Code is released under the MIT License. The data are U.S. government public data.
