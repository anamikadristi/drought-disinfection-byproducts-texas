"""Pilot: does drought raise disinfection byproduct (DBP) violations in Texas drinking water?

Inputs (data/):
  tx_cws_systems_2015_2025.csv         Texas community water systems (EPA SDWIS via the drinking-water-violations-tx-la project)
  tx_dbp_violations_2015_2025.csv      Health-based DBP violations (TTHM, HAA5 MCL; MRDL), Texas, 2015-2025
  usdm_tx_county_weekly_2014_2025.csv  U.S. Drought Monitor, percent of county area in each category (weekly)
  nclimdiv_tx_county_tmean_1991_2025.csv NOAA nClimDiv county monthly mean temperature (deg F)

Design: system-quarter panel, 2016Q1-2025Q4. Outcome = number of TTHM or HAA5 MCL violations whose
compliance period begins in that quarter. Exposure = county percent area in drought (D1 or worse),
averaged over the preceding four quarters (the window of the Stage 2 locational running annual average).
Models: Poisson GLM with year and season terms, system size and source-water type, county temperature
anomaly, and county-clustered standard errors; run for all systems and by source-water type.
"""
import os, numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D, T, F = P + "/data/", P + "/outputs/tables/", P + "/outputs/figures/"

# ---------------- drought: county-quarter mean % area D1+, then trailing 4-quarter mean
u = pd.read_csv(D + "usdm_tx_county_weekly_2014_2025.csv", dtype={"FIPS": str})
u["date"] = pd.to_datetime(u.MapDate.astype(str)); u = u[u.date >= "2014-01-01"]
u["q"] = u.date.dt.to_period("Q")
dq = u.groupby(["FIPS", "q"]).D1.mean().rename("d1").reset_index().sort_values(["FIPS", "q"])
dq["d1_4q"] = dq.groupby("FIPS").d1.transform(lambda x: x.rolling(4, min_periods=4).mean())
dq["d1_prev4q"] = dq.groupby("FIPS").d1_4q.shift(4)

# ---------------- temperature: county-quarter anomaly vs 1991-2020 same-quarter mean
t = pd.read_csv(D + "nclimdiv_tx_county_tmean_1991_2025.csv", dtype={"cty3": str})
tl = t.melt(id_vars=["cty3", "year"], var_name="m", value_name="tf"); tl = tl[tl.tf > -90]
tl["month"] = tl.m.str[1:].astype(int); tl["qn"] = (tl.month - 1) // 3 + 1
tq = tl.groupby(["cty3", "year", "qn"]).tf.mean().reset_index()
clim = tq[(tq.year >= 1991) & (tq.year <= 2020)].groupby(["cty3", "qn"]).tf.mean().rename("tclim").reset_index()
tq = tq.merge(clim, on=["cty3", "qn"]); tq["tanom"] = tq.tf - tq.tclim
tq["FIPS"] = "48" + tq.cty3; tq["q"] = pd.PeriodIndex.from_fields(year=tq.year, quarter=tq.qn, freq="Q")

# ---------------- system-quarter panel
s = pd.read_csv(D + "tx_cws_systems_2015_2025.csv", dtype={"county_fips": str})
s = s.dropna(subset=["county_fips"]); s["county_fips"] = s.county_fips.str.zfill(5)
s["src"] = s.source.map({"Surface water": "Surface water", "Surface water (purchased)": "Surface water",
                         "Groundwater": "Groundwater", "Groundwater (purchased)": "Groundwater"})
s["own_surface"] = (s.source == "Surface water").astype(int)
qs = pd.period_range("2016Q1", "2025Q4", freq="Q")
rows = []
for r in s.itertuples():
    st, en = pd.Period(pd.Timestamp(r.start), "Q"), pd.Period(pd.Timestamp(r.end), "Q")
    for q in qs:
        if st <= q <= en: rows.append((r.PWSID, r.county_fips, r.src, r.source, r.size_class, q))
pan = pd.DataFrame(rows, columns=["PWSID", "FIPS", "src", "source", "size_class", "q"])
v = pd.read_csv(D + "tx_dbp_violations_2015_2025.csv", dtype=str)
v = v[(v.VIOLATION_CATEGORY_CODE == "MCL") & v.contaminant.isin(["TTHM", "Total Haloacetic Acids (HAA5)"])]
v["q"] = pd.to_datetime(v.begin).dt.to_period("Q")
vc = v.groupby(["PWSID", "q"]).size().rename("n_viol").reset_index()
pan = pan.merge(vc, on=["PWSID", "q"], how="left").fillna({"n_viol": 0})
pan = pan.merge(dq, on=["FIPS", "q"], how="left").merge(tq[["FIPS", "q", "tanom"]], on=["FIPS", "q"], how="left")
pan["year"] = pan.q.dt.year; pan["season"] = pan.q.dt.quarter
pan["drought10"] = pan.d1_4q / 10; pan["prev10"] = pan.d1_prev4q / 10
pan = pan.dropna(subset=["drought10", "tanom"])
pan["any_viol"] = (pan.n_viol > 0).astype(int)

# ---------------- descriptive tables
desc = pan.groupby("src").agg(systems=("PWSID", "nunique"), system_quarters=("PWSID", "size"),
                              violations=("n_viol", "sum"), mean_drought_pct=("d1_4q", "mean")).reset_index()
desc["viol_per_100_sys_yr"] = 100 * desc.violations / (desc.system_quarters / 4)
desc.to_csv(T + "t1_descriptive_by_source.csv", index=False)
pan["drought_tertile"] = pd.qcut(pan.d1_4q.rank(method="first"), 3, labels=["Low", "Middle", "High"])
ter = pan.groupby(["src", "drought_tertile"], observed=True).agg(system_quarters=("PWSID", "size"), violations=("n_viol", "sum"),
                                                                 drought_pct=("d1_4q", "mean")).reset_index()
ter["viol_per_100_sys_yr"] = 100 * ter.violations / (ter.system_quarters / 4)
ter.to_csv(T + "t2_rates_by_drought_tertile.csv", index=False)

# ---------------- models
def fit(df, label, formula):
    cols = ["n_viol", "drought10", "tanom", "year", "season", "size_class", "source", "FIPS"] + (["prev10"] if "prev10" in formula else [])
    df = df.dropna(subset=cols)
    m = smf.glm(formula, data=df, family=sm.families.Poisson()).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(df.FIPS)[0]})
    out = []
    for term in ["drought10", "prev10", "tanom"]:
        if term in m.params:
            b, se = m.params[term], m.bse[term]
            out.append(dict(model=label, term=term, IRR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=m.pvalues[term], n=int(m.nobs)))
    return out
base = "n_viol ~ drought10 + tanom + C(year) + C(season) + C(size_class)"
res = []
res += fit(pan, "All systems", base + " + C(source)")
res += fit(pan[pan.src == "Surface water"], "Surface-water systems", base + " + C(source)")
res += fit(pan[pan.src == "Groundwater"], "Groundwater systems", base + " + C(source)")
res += fit(pan.dropna(subset=["prev10"]), "All systems, current vs previous year drought", base + " + prev10 + C(source)")
# county fixed effects (within-county change over time)
cfe = pan.groupby("FIPS").n_viol.transform("sum") > 0
res += fit(pan[cfe], "All systems, county fixed effects", base + " + C(source) + C(FIPS)")
res = pd.DataFrame(res); res.to_csv(T + "t3_poisson_irr.csv", index=False)

# ---------------- statewide quarterly series for the figure
sq = pan[pan.src == "Surface water"].groupby("q").agg(viol=("n_viol", "sum"), sys=("PWSID", "nunique")).reset_index()
sq["rate"] = 100 * sq.viol / sq.sys
st = pan.groupby("q").d1_4q.mean().reset_index()
sq = sq.merge(st, on="q"); sq.to_csv(T + "t4_statewide_quarterly.csv", index=False)

fig, ax = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw={"width_ratios": [1.6, 1]})
x = sq.q.dt.to_timestamp()
ax[0].fill_between(x, sq.d1_4q, color="#E8B04A", alpha=0.45, label="Drought (D1+), % county area, trailing 4 quarters")
ax[0].set_ylabel("Drought area (%)"); ax2 = ax[0].twinx()
ax2.plot(x, sq.rate, color="#1F3A5F", lw=2, label="DBP MCL violations per 100 surface-water systems")
ax2.set_ylabel("Violations per 100 systems per quarter"); ax[0].set_title("Texas, 2016-2025")
h1, l1 = ax[0].get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax[0].legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right")
r = res[res.term == "drought10"].iloc[::-1]
ax[1].errorbar(r.IRR, range(len(r)), xerr=[r.IRR - r.lo, r.hi - r.IRR], fmt="o", color="#1F3A5F", capsize=3)
ax[1].axvline(1, color="grey", ls="--", lw=1); ax[1].set_yticks(range(len(r))); ax[1].set_yticklabels(r.model, fontsize=7)
ax[1].set_xlabel("Rate ratio per +10 pts drought area"); ax[1].set_title("Poisson models (95% CI)")
plt.tight_layout(); plt.savefig(F + "drought_dbp_texas.png", dpi=200)
print(desc.to_string(index=False)); print(ter.to_string(index=False)); print(res.round(3).to_string(index=False))
