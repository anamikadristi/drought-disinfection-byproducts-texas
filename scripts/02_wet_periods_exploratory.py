"""Exploratory follow-up: are DBP violations higher after unusually WET periods (organic carbon flushing),
and after wet periods that follow drought (drought-to-flood whiplash)?
Exposure: county precipitation anomaly (z-score vs 1991-2020, same quarter) in the current and previous quarter.
Same panel, covariates and Poisson models as 01_analysis.py. Labeled exploratory: tested after the drought result."""
import os, runpy, numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
g = runpy.run_path(P + "/scripts/01_analysis.py")          # rebuilds panel and saves the drought outputs
pan = g["pan"]; T = P + "/outputs/tables/"
pr = pd.read_csv(P + "/data/nclimdiv_tx_county_precip_1991_2025.csv", dtype={"cty3": str})
pl = pr.melt(id_vars=["cty3", "year"], var_name="m", value_name="pin"); pl = pl[pl.pin > -9]
pl["qn"] = (pl.m.str[1:].astype(int) - 1) // 3 + 1
pq = pl.groupby(["cty3", "year", "qn"]).pin.sum().reset_index()
c = pq[(pq.year >= 1991) & (pq.year <= 2020)].groupby(["cty3", "qn"]).pin.agg(["mean", "std"]).reset_index()
pq = pq.merge(c, on=["cty3", "qn"]); pq["pz"] = (pq.pin - pq["mean"]) / pq["std"]
pq["FIPS"] = "48" + pq.cty3; pq["q"] = pd.PeriodIndex.from_fields(year=pq.year, quarter=pq.qn, freq="Q")
pq = pq.sort_values(["FIPS", "q"]); pq["pz_prev"] = pq.groupby("FIPS").pz.shift(1)
pan = pan.merge(pq[["FIPS", "q", "pz", "pz_prev"]], on=["FIPS", "q"], how="left").dropna(subset=["pz", "pz_prev"])
pan["wet_prev"] = (pan.pz_prev > 1).astype(int)
pan["whiplash"] = ((pan.pz_prev > 1) & (pan.d1_prev4q > 30)).astype(int)
def fit(df, label, f, terms):
    df = df.dropna(subset=["size_class", "source", "tanom"])
    m = smf.glm(f, data=df, family=sm.families.Poisson()).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(df.FIPS)[0]})
    return [dict(model=label, term=t, IRR=np.exp(m.params[t]), lo=np.exp(m.params[t] - 1.96 * m.bse[t]), hi=np.exp(m.params[t] + 1.96 * m.bse[t]), p=m.pvalues[t], n=int(m.nobs)) for t in terms]
b = "n_viol ~ tanom + C(year) + C(season) + C(size_class) + C(source)"
out = []
for lab, df in [("All systems", pan), ("Surface-water systems", pan[pan.src == "Surface water"]), ("Groundwater systems", pan[pan.src == "Groundwater"])]:
    out += fit(df, lab + ": precipitation z (current, previous quarter)", b + " + pz + pz_prev", ["pz", "pz_prev"])
    out += fit(df, lab + ": very wet previous quarter (z>1)", b + " + wet_prev", ["wet_prev"])
    out += fit(df.dropna(subset=["d1_prev4q"]), lab + ": wet quarter after drought year (whiplash)", b + " + wet_prev + whiplash + drought10", ["wet_prev", "whiplash"])
out = pd.DataFrame(out); out.to_csv(T + "t5_wet_periods_exploratory.csv", index=False)
print(out.round(3).to_string(index=False))
print("whiplash system-quarters:", int(pan.whiplash.sum()), "| wet_prev:", int(pan.wet_prev.sum()))
