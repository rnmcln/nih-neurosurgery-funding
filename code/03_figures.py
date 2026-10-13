#!/usr/bin/env python3
"""
03_figures.py
=============
Publication figures with a single cohesive visual identity:
  * one typeface (DejaVu Sans) and one colour palette across every figure;
  * no in-figure titles -- all description is carried by the figure legends;
  * keys placed beneath the plotting area so they never overlap data;
  * a clean, minimal R/ggplot-style theme (white panel, light horizontal
    gridlines, no top/right spine).

Reads the locally generated outputs/ directory when present, otherwise the
bundled aggregate results in data/aggregate/. Writes PNG (300 dpi) and TIFF.
"""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCAL_OUT = os.path.join(ROOT, "outputs")
BUNDLED_OUT = os.path.join(ROOT, "data", "aggregate")
DEFAULT_OUT = (LOCAL_OUT
               if os.path.exists(os.path.join(LOCAL_OUT, "analysis_output.json"))
               else BUNDLED_OUT)
OUT = os.environ.get("REPORTER_OUTPUT_DIR", DEFAULT_OUT)
FIGDIR = os.environ.get("REPORTER_FIGURE_DIR", os.path.join(ROOT, "figures"))
os.makedirs(FIGDIR, exist_ok=True)
A = json.load(open(os.path.join(OUT, "analysis_output.json")))
M = json.load(open(os.path.join(OUT, "pubmed_mesh_lag.json")))

# ---------------------------------------------------------------- palette
INK = "#33363B"        # charcoal: text, axes
GRID = "#E6E6E6"       # light horizontal gridlines
PRIMARY = "#33628F"    # steel blue
PRIMARY_L = "#9DB6CC"  # light steel blue
SECOND = "#E0A22A"     # goldenrod
THIRD = "#4E937A"      # sage green
FOURTH = "#B5563F"     # terracotta
GREY = "#9AA0A6"       # cool grey

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5, "text.color": INK,
    "axes.labelsize": 8.5, "axes.labelcolor": INK, "axes.edgecolor": "#BFBFBF",
    "axes.linewidth": 0.8, "axes.titlesize": 8.5,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 7.6, "figure.dpi": 300,
})
YEARS = list(range(2015, 2026))


def theme(ax, grid_axis="y"):
    """Apply the shared minimal R-style theme to an axis."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#BFBFBF")
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, color="#BFBFBF")


def legend_below(ax, ncol=3, y=-0.30):
    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol,
                    frameon=False, handlelength=1.6, columnspacing=1.6,
                    handletextpad=0.5)
    return leg


def tag(ax, s):
    ax.text(-0.13, 1.04, s, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", ha="left", color=INK)


def save(fig, name):
    for ext in ("png", "tiff"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"), dpi=300,
                    bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"} if ext == "tiff" else {})
    plt.close(fig)
    print(f"  saved {name}")


# =========================================================================
# FIGURE 1 -- annual funding (single panel; budget share moved to SDC)
# =========================================================================
nominal = [A["fy_nominal"][str(y)] / 1e6 for y in YEARS]
real = [A["fy_real"][str(y)] / 1e6 for y in YEARS]
fit_yrs = list(range(2015, 2024))
logf = np.log([A["fy_real"][str(y)] for y in fit_yrs])
sl, ic, r, p, se = stats.linregress(fit_yrs, logf)
resid = logf - (ic + sl * np.array(fit_yrs))
s_err = np.sqrt(np.sum(resid ** 2) / (len(fit_yrs) - 2))
xb = np.mean(fit_yrs)
sxx = np.sum((np.array(fit_yrs) - xb) ** 2)
px = np.array(range(2015, 2026))
proj = np.exp(ic + sl * px) / 1e6
piw = np.array([stats.t.ppf(0.975, len(fit_yrs) - 2) * s_err *
                math.sqrt(1 + 1 / len(fit_yrs) + (x - xb) ** 2 / sxx) for x in px])
plo = np.exp((ic + sl * px) - piw) / 1e6
phi = np.exp((ic + sl * px) + piw) / 1e6

fig, ax = plt.subplots(figsize=(5.4, 3.3))
x = np.arange(len(YEARS))
w = 0.40
ax.fill_between(x, plo, phi, color=GREY, alpha=0.18, zorder=1,
                label="FY2015–FY2023 trend (95% prediction interval)")
ax.plot(x, proj, color=GREY, lw=1.1, ls=(0, (5, 2)), zorder=2)
ax.bar(x - w / 2, nominal, w, color=PRIMARY_L, zorder=3, label="Nominal dollars")
ax.bar(x + w / 2, real, w, color=PRIMARY, zorder=3,
       label="Constant 2025 dollars (CPI-U)")
theme(ax)
ax.set_xticks(x)
ax.set_xticklabels([f"FY{str(y)[2:]}" for y in YEARS])
ax.set_ylabel("Annual NIH funding (US$ millions)")
ax.set_ylim(0, 360)
legend_below(ax, ncol=2, y=-0.26)
save(fig, "figure1_funding")

# =========================================================================
# FIGURE 2 -- institutional concentration (2 panels)
# =========================================================================
con = A["concentration"]
gini = [con[str(y)]["gini"] for y in YEARS]
top5 = [con[str(y)]["top5_pct"] for y in YEARS]
top10 = [con[str(y)]["top10_pct"] for y in YEARS]
ndep = [con[str(y)]["n_departments"] for y in YEARS]
ret = [A["sustainability"]["yoy_retention"][str(y)] for y in range(2015, 2025)]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.3, 3.5))
xl = [f"FY{str(y)[2:]}" for y in YEARS]

ax1.plot(x, top10, marker="o", ms=4, color=SECOND, lw=1.6, label="Top-10 share")
ax1.plot(x, top5, marker="s", ms=3.6, color=FOURTH, lw=1.6, label="Top-5 share")
ax1.set_ylabel("Share of annual funding (%)")
ax1.set_ylim(0, 80)
theme(ax1)
ax1b = ax1.twinx()
ax1b.plot(x, gini, marker="D", ms=3.6, color=PRIMARY, lw=1.8, label="Gini coefficient")
ax1b.set_ylabel("Gini coefficient", color=PRIMARY)
ax1b.set_ylim(0, 1.0)
ax1b.tick_params(axis="y", colors=PRIMARY)
ax1b.spines["top"].set_visible(False)
ax1b.spines["right"].set_color(PRIMARY)
ax1.set_xticks(x); ax1.set_xticklabels(xl, rotation=45, ha="right")
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax1b.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, loc="upper center", bbox_to_anchor=(0.5, -0.34),
           ncol=3, frameon=False, handlelength=1.6, columnspacing=1.3)
tag(ax1, "A")

ax2.bar(x, ndep, color=PRIMARY_L, width=0.66, label="Funded departments")
for xi, n in zip(x, ndep):
    ax2.text(xi, n + 1.3, str(n), ha="center", va="bottom", fontsize=6.6, color=INK)
ax2.set_ylabel("Funded departments (n)")
ax2.set_ylim(0, 72)
theme(ax2)
ax2b = ax2.twinx()
ax2b.plot(x[:-1], ret, marker="o", ms=3.6, color=THIRD, lw=1.7,
          label="Retention to next fiscal year")
ax2b.set_ylabel("Retention to next FY (%)", color=THIRD)
ax2b.set_ylim(0, 100)
ax2b.tick_params(axis="y", colors=THIRD)
ax2b.spines["top"].set_visible(False)
ax2b.spines["right"].set_color(THIRD)
ax2.set_xticks(x); ax2.set_xticklabels(xl, rotation=45, ha="right")
h1, l1 = ax2.get_legend_handles_labels()
h2, l2 = ax2b.get_legend_handles_labels()
ax2.legend(h1 + h2, l1 + l2, loc="upper center", bbox_to_anchor=(0.5, -0.34),
           ncol=2, frameon=False, handlelength=1.6, columnspacing=1.3)
tag(ax2, "B")
fig.subplots_adjust(wspace=0.42)
save(fig, "figure2_concentration")

# =========================================================================
# FIGURE 3 -- keyword-attributed disease-domain funding
# =========================================================================
sA = A["subspecialty_real_funding"]
sB = A["subspecialty_real_funding_B"]
order = ["Neuro-oncology", "Spine", "Cerebrovascular",
         "Functional/epilepsy", "Trauma/TBI", "Hydrocephalus",
         "Peripheral nerve"]
vA = [sA[s] / 1e6 for s in order]
vB = [sB[s] / 1e6 for s in order]

fig, ax = plt.subplots(figsize=(5.2, 3.8))
y = np.arange(len(order))
hh = 0.38
ax.barh(y + hh / 2, vA, hh, color=PRIMARY, zorder=3,
        label="Primary classifier (expanded text)")
ax.barh(y - hh / 2, vB, hh, color=PRIMARY_L, zorder=3,
        label="Robustness classifier (title + spending category)")
for yi, v in zip(y, vA):
    ax.text(v + 9, yi + hh / 2, f"{v:,.0f}", va="center", fontsize=6.8, color=INK)
ax.set_yticks(y)
ax.set_yticklabels(order)
ax.invert_yaxis()
ax.set_xlabel("Cumulative real funding, FY2015–FY2025 (US$ millions)")
ax.set_xlim(0, max(vA) * 1.16)
theme(ax, grid_axis="x")
legend_below(ax, ncol=1, y=-0.26)
save(fig, "figure3_domain")

# =========================================================================
# FIGURE 4 -- training infrastructure (2 panels)
# =========================================================================
pp = A["pi_pipeline"]
tm = pp["t32_adjusted_model"]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.3, 3.4))

# 4A -- transition by entry mechanism
groups = ["Fellowship\n(F30/F31/F32)", "Mentored career\ndevelopment (K)"]
rates = [pp["fellowship_rate"], pp["mentored_K_rate"]]
ns = [pp["fellowship_cohort"], pp["mentored_K_cohort"]]
trans = [pp["fellowship_transitioned"], pp["mentored_K_transitioned"]]
ax1.bar([0, 1], rates, 0.56, color=[PRIMARY_L, PRIMARY], zorder=3)
for i, (rt, n, tr) in enumerate(zip(rates, ns, trans)):
    ax1.text(i, rt + 1.6, f"{rt:.0f}%", ha="center", va="bottom",
             fontsize=8.5, fontweight="bold", color=INK)
    ax1.text(i, rt + 6.2, f"{tr}/{n}", ha="center", va="bottom",
             fontsize=6.8, color=GREY)
ax1.set_xticks([0, 1]); ax1.set_xticklabels(groups)
ax1.set_ylabel("Later held an independent R-type award (%)")
ax1.set_ylim(0, 52)
theme(ax1)
tag(ax1, "A")

# 4B -- T32 forest plot: unadjusted vs adjusted incidence rate ratio.
# Short y-tick labels keep all text clear of the data; the covariate set
# for the adjusted model is given in the figure legend.
labels = ["Unadjusted", "Adjusted"]
irrs = [tm["unadjusted_irr"], tm["adjusted_irr"]]
los = [tm["unadjusted_irr_ci"][0], tm["adjusted_irr_ci"][0]]
his = [tm["unadjusted_irr_ci"][1], tm["adjusted_irr_ci"][1]]
yy = [1, 0]
cols = [FOURTH, PRIMARY]
ax2.axvline(1.0, color=GREY, lw=1.0, ls=(0, (4, 2)), zorder=1)
for yi, irr, lo, hi, c in zip(yy, irrs, los, his, cols):
    ax2.plot([lo, hi], [yi, yi], color=c, lw=2.0, zorder=3)
    ax2.plot(irr, yi, "o", color=c, ms=8, zorder=4)
    ax2.text(irr, yi + 0.20, f"IRR {irr:.2f}", ha="center", fontsize=7.4,
             color=INK)
ax2.set_xscale("log")
ax2.set_xlim(0.3, 9)
ax2.set_xticks([0.5, 1, 2, 4, 8])
ax2.set_xticklabels(["0.5", "1", "2", "4", "8"])
ax2.set_yticks(yy)
ax2.set_yticklabels(labels)
ax2.set_ylim(-0.5, 1.5)
ax2.set_xlabel("Incidence rate ratio (distinct R01 PIs)")
theme(ax2, grid_axis="x")
tag(ax2, "B")
fig.subplots_adjust(wspace=0.55)
save(fig, "figure4_training")

# =========================================================================
# FIGURE 5 -- MeSH retrieval coverage
# =========================================================================
yrs = list(range(2015, 2026))
gen = [M["years"][str(y)]["mesh_to_keyword_ratio"] for y in yrs]
jrn = [M["years"][str(y)]["journal_mesh_ratio"] for y in yrs]
xi = np.arange(len(yrs))
fig, ax = plt.subplots(figsize=(5.4, 3.4))
ax.axvspan(7.5, 10.5, color=GREY, alpha=0.10, zorder=0)
ax.plot(xi[:9], gen[:9], marker="o", ms=4.4, color=PRIMARY, lw=1.9,
        label="General PubMed (MeSH ÷ keyword-augmented)")
ax.plot(xi[8:], gen[8:], marker="o", ms=3.6, color=PRIMARY, lw=1.1,
        ls=(0, (2, 1.5)), alpha=0.7)
ax.plot(xi[:9], jrn[:9], marker="s", ms=4.0, color=SECOND, lw=1.9,
        label="Nine-journal subset (MeSH-indexed ÷ all)")
ax.plot(xi[8:], jrn[8:], marker="s", ms=3.4, color=SECOND, lw=1.1,
        ls=(0, (2, 1.5)), alpha=0.7)
ax.text(9.0, 0.035, "incomplete\nindexing", ha="center", fontsize=6.4,
        color=GREY, style="italic")
ax.set_xticks(xi)
ax.set_xticklabels([f"'{str(y)[2:]}" for y in yrs])
ax.set_xlabel("Publication year")
ax.set_ylabel("MeSH retrieval coverage ratio")
ax.set_ylim(0, 0.75)
theme(ax)
legend_below(ax, ncol=1, y=-0.26)
save(fig, "figure5_mesh")

# =========================================================================
# SDC FIGURE 1 -- Lorenz curve, FY2025
# =========================================================================
import pandas as pd
cohort = pd.read_csv(os.path.join(OUT, "deduplicated_cohort.csv"))
fy25 = cohort[cohort["Fiscal Year"] == 2025]
inst = np.sort(fy25.groupby("Organization Name")["Total Cost"].sum().dropna().values)
inst = inst[inst > 0]
cum = np.insert(np.cumsum(inst) / inst.sum(), 0, 0)
xx = np.linspace(0, 1, len(cum))
fig, ax = plt.subplots(figsize=(3.9, 3.9))
ax.plot([0, 1], [0, 1], color=GREY, lw=1.0, ls=(0, (4, 2)),
        label="Line of equality")
ax.plot(xx, cum, color=PRIMARY, lw=2.0, label="FY2025 funding")
ax.fill_between(xx, cum, xx, color=PRIMARY, alpha=0.14)
ax.set_xlabel("Cumulative share of funded departments")
ax.set_ylabel("Cumulative share of funding")
ax.text(0.06, 0.86, f"Gini = {con['2025']['gini']:.3f}", fontsize=8.5,
        fontweight="bold", color=INK)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
theme(ax, grid_axis="both")
legend_below(ax, ncol=2, y=-0.20)
save(fig, "sfigure1_lorenz")

# =========================================================================
# SDC FIGURE 2 -- NIH-wide benchmark, indexed real funding
# =========================================================================
dept = np.array([A["fy_real"][str(y)] for y in YEARS], dtype=float)
rpg = np.array([A["nih_rpg_real"][str(y)] for y in YEARS], dtype=float)
app = np.array([A["nih_appropriation_real"][str(y)] for y in YEARS], dtype=float)
dept_i, rpg_i, app_i = (100 * dept / dept[0], 100 * rpg / rpg[0],
                        100 * app / app[0])
fig, ax = plt.subplots(figsize=(5.4, 3.2))
ax.plot(x, dept_i, marker="o", ms=4.2, color=PRIMARY, lw=1.9,
        label="Neurosurgery-department cohort")
ax.plot(x, rpg_i, marker="s", ms=3.8, color=SECOND, lw=1.7,
        label="NIH-wide research project grants")
ax.plot(x, app_i, marker="D", ms=3.5, color=THIRD, lw=1.7,
        label="Total NIH appropriation")
ax.axhline(100, color=GREY, lw=0.8, ls=(0, (3, 2)))
ax.set_xticks(x)
ax.set_xticklabels([f"FY{str(y)[2:]}" for y in YEARS])
ax.set_ylabel("Real funding index (FY2015 = 100)")
ax.set_ylim(80, 220)
theme(ax)
legend_below(ax, ncol=1, y=-0.28)
save(fig, "sfigure2_benchmark")

# =========================================================================
# SDC FIGURE 3 -- deduplication sensitivity (raw vs deduplicated real funding)
# =========================================================================
ds = A["dedup_sensitivity"]
raw_real = [ds["raw_real"][str(y)] / 1e6 for y in YEARS]
ddup_real = [A["fy_real"][str(y)] / 1e6 for y in YEARS]
fig, ax = plt.subplots(figsize=(5.4, 3.2))
ax.plot(x, raw_real, marker="o", ms=4.0, color=GREY, lw=1.7,
        label="Before deduplication (8,672 export rows)")
ax.plot(x, ddup_real, marker="D", ms=3.8, color=PRIMARY, lw=1.9,
        label="After deduplication (8,480 records)")
ax.set_xticks(x)
ax.set_xticklabels([f"FY{str(y)[2:]}" for y in YEARS])
ax.set_ylabel("Annual funding, constant 2025 US$ (millions)")
ax.set_ylim(0, 320)
theme(ax)
legend_below(ax, ncol=1, y=-0.26)
save(fig, "sfigure3_dedup")

print("\nAll figures written to", FIGDIR)
