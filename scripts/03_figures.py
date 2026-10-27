"""Summary figure: drought vs wet-period rate ratios for DBP violations (Texas community water systems, 2016-2025)."""
import os, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); T = P + "/outputs/tables/"
d = pd.read_csv(T + "t3_poisson_irr.csv"); w = pd.read_csv(T + "t5_wet_periods_exploratory.csv")
rows = []
for src in ["Surface-water systems", "Groundwater systems", "All systems"]:
    a = d[(d.model == src) & (d.term == "drought10")].iloc[0]; rows.append((src, "Drought, +10 pts county area (past year)", a))
    b = w[(w.model == f"{src}: very wet previous quarter (z>1)")].iloc[0]; rows.append((src, "Very wet previous quarter (rain z > 1)", b))
    c = w[(w.model == f"{src}: wet quarter after drought year (whiplash)") & (w.term == "whiplash")].iloc[0]; rows.append((src, "Extra effect: wet quarter after a drought year", c))
fig, ax = plt.subplots(figsize=(8.2, 5.6)); y = 0; ticks = []; labs = []
cols = {"Drought, +10 pts county area (past year)": "#E8B04A", "Very wet previous quarter (rain z > 1)": "#3A7CA5", "Extra effect: wet quarter after a drought year": "#1F3A5F"}
for i, (src, lab, r) in enumerate(rows):
    if i and i % 3 == 0: y += 0.8
    ax.errorbar(r.IRR, y, xerr=[[r.IRR - r.lo], [r.hi - r.IRR]], fmt="o", color=cols[lab], capsize=3, label=lab if i < 3 else None)
    ticks.append(y); labs.append(src if i % 3 == 1 else ""); y += 1
ax.axvline(1, color="grey", ls="--", lw=1); ax.set_yticks(ticks); ax.set_yticklabels(labs, fontsize=9); ax.invert_yaxis()
ax.set_xlabel("Rate ratio for DBP (TTHM/HAA5) MCL violations, 95% CI"); ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1, frameon=False)
ax.set_title("Texas community water systems, 2016-2025: wet periods, not drought, track DBP violations", fontsize=10)
plt.tight_layout(); plt.savefig(P + "/outputs/figures/dbp_drought_vs_wet_texas.png", dpi=200); print("ok")
