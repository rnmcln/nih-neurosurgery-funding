#!/usr/bin/env python3
"""Generate publication figures with a restrained, high-data-ink visual system."""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCAL_OUT = os.path.join(ROOT, "analysis", "outputs")
BUNDLED_OUT = os.path.join(ROOT, "data", "aggregate")
DEFAULT_OUT = (LOCAL_OUT
               if os.path.exists(os.path.join(LOCAL_OUT, "analysis_output.json"))
               else BUNDLED_OUT)
OUT = os.environ.get("REPORTER_OUTPUT_DIR", DEFAULT_OUT)
FIGDIR = os.environ.get("REPORTER_FIGURE_DIR", os.path.join(ROOT, "figures"))
os.makedirs(FIGDIR, exist_ok=True)
with open(os.path.join(OUT, "analysis_output.json"), encoding="utf-8") as f:
    A = json.load(f)
with open(os.path.join(OUT, "pubmed_mesh_lag.json"), encoding="utf-8") as f:
    M = json.load(f)

# Quiet, colour-blind-safe palette. Colour identifies data, never decoration.
INK = "#24272B"
MUTED = "#6F767D"
PALE = "#AEB4BA"
RULE = "#D9DCDF"
PRIMARY = "#315F86"
PRIMARY_L = "#8EA8BC"
SECOND = "#B07D1D"
THIRD = "#3D7C69"
FOURTH = "#A84F3A"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.3, "text.color": INK,
    "axes.labelsize": 8.3, "axes.labelcolor": INK,
    "axes.edgecolor": PALE, "axes.linewidth": 0.65,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.labelsize": 7.6, "ytick.labelsize": 7.6,
    "legend.fontsize": 7.2, "figure.dpi": 300, "savefig.facecolor": "white",
})

YEARS = list(range(2015, 2026))
X = np.arange(len(YEARS))
FY_LABELS = [f"FY{str(y)[2:]}" for y in YEARS]


def theme(ax, grid_axis="y"):
    """Minimal axis treatment: data first, scaffolding second."""
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(PALE)
    ax.spines["bottom"].set_linewidth(0.65)
    if grid_axis:
        ax.grid(axis=grid_axis, color=RULE, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.tick_params(axis="x", length=2.5, width=0.6, color=PALE, pad=4)


def tag(ax, label):
    ax.text(-0.11, 1.035, label, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left", color=INK)


def sparse_year_ticks(ax, top_row=False):
    ticks = [0, 2, 4, 6, 8, 10]
    ax.set_xticks(ticks)
    ax.set_xticklabels([] if top_row else [FY_LABELS[i] for i in ticks])


def save(fig, name):
    for ext in ("png", "tiff"):
        kwargs = {"compression": "tiff_lzw"} if ext == "tiff" else {}
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"), dpi=300,
                    bbox_inches="tight", pad_inches=0.08, pil_kwargs=kwargs)
    plt.close(fig)
    print(f"  saved {name}")


# FIGURE 1 -- annual funding and pre-2024 trend
nominal = np.array([A["fy_nominal"][str(y)] / 1e6 for y in YEARS])
real = np.array([A["fy_real"][str(y)] / 1e6 for y in YEARS])
fit_yrs = np.array(range(2015, 2024))
logf = np.log([A["fy_real"][str(y)] for y in fit_yrs])
sl, ic, r, p, se = stats.linregress(fit_yrs, logf)
resid = logf - (ic + sl * fit_yrs)
s_err = np.sqrt(np.sum(resid ** 2) / (len(fit_yrs) - 2))
xb = np.mean(fit_yrs)
sxx = np.sum((fit_yrs - xb) ** 2)
px = np.array(YEARS)
proj = np.exp(ic + sl * px) / 1e6
piw = np.array([stats.t.ppf(0.975, len(fit_yrs) - 2) * s_err
                * math.sqrt(1 + 1 / len(fit_yrs) + (year - xb) ** 2 / sxx)
                for year in px])
plo = np.exp((ic + sl * px) - piw) / 1e6
phi = np.exp((ic + sl * px) + piw) / 1e6

fig, ax = plt.subplots(figsize=(6.15, 3.25))
ax.fill_between(X, plo, phi, color=PALE, alpha=0.18, linewidth=0, zorder=1)
ax.plot(X, proj, color=MUTED, lw=1.0, ls=(0, (4, 2)), zorder=2)
nominal_line, = ax.plot(X, nominal, color=PRIMARY_L, lw=1.35, marker="o",
                        ms=3.0, mfc="white", mew=0.8, zorder=3)
real_line, = ax.plot(X, real, color=PRIMARY, lw=2.0, marker="o", ms=3.2,
                     zorder=4)
legend = ax.legend(
    [real_line, nominal_line], ["Constant 2025 dollars", "Nominal dollars"],
    loc="upper left", bbox_to_anchor=(0.01, 0.995), borderaxespad=0,
    frameon=False, handlelength=2.5, handletextpad=0.7, labelspacing=0.55,
    fontsize=7.5,
)
for label, color in zip(legend.get_texts(), [PRIMARY, PRIMARY_L]):
    label.set_color(color)
ax.text(6.85, 352, "FY2015–FY2023 trend\n95% prediction interval",
        color=MUTED, fontsize=6.9, ha="left", va="top")
ax.text(10.0, real[-1] - 14, f"FY2025  ${real[-1]:.0f}M", color=INK,
        fontsize=7.1, ha="right", va="top")
ax.set(xlim=(-0.35, 10.35), ylim=(80, 365),
       ylabel="Annual NIH funding (US$ millions)")
ax.set_xticks(X)
ax.set_xticklabels(FY_LABELS)
theme(ax)
fig.subplots_adjust(left=0.10, right=0.985, bottom=0.16, top=0.97)
save(fig, "figure1_funding")


# FIGURE 2 -- institutional concentration as four aligned small multiples
con = A["concentration"]
gini = np.array([con[str(y)]["gini"] for y in YEARS])
top5 = np.array([con[str(y)]["top5_pct"] for y in YEARS])
top10 = np.array([con[str(y)]["top10_pct"] for y in YEARS])
ndep = np.array([con[str(y)]["n_departments"] for y in YEARS])
ret = np.array([A["sustainability"]["yoy_retention"][str(y)]
                for y in range(2015, 2025)])
fig, axs = plt.subplots(2, 2, figsize=(7.35, 4.65), sharex=True)
ax1, ax2, ax3, ax4 = axs.ravel()

ax1.plot(X, top10, color=PRIMARY, lw=1.75)
ax1.plot(X, top5, color=PRIMARY_L, lw=1.45)
ax1.scatter(X, top10, s=10, color=PRIMARY, zorder=3)
ax1.scatter(X, top5, s=9, color=PRIMARY_L, zorder=3)
ax1.text(10.18, top10[-1], f"Top 10  {top10[-1]:.1f}%", color=PRIMARY,
         fontsize=6.9, va="center")
ax1.text(10.18, top5[-1], f"Top 5  {top5[-1]:.1f}%", color=PRIMARY_L,
         fontsize=6.9, va="center")
ax1.set(xlim=(-0.25, 11.1), ylim=(30, 65), ylabel="Share of annual funding (%)")
theme(ax1); sparse_year_ticks(ax1, top_row=True); tag(ax1, "A")

ax2.plot(X, gini, color=INK, lw=1.55)
ax2.scatter(X, gini, s=10, color=INK, zorder=3)
ax2.text(10.18, gini[-1], f"{gini[-1]:.3f}", color=INK,
         fontsize=6.9, va="center")
ax2.set(xlim=(-0.25, 11.1), ylim=(0.52, 0.59), ylabel="Gini coefficient")
ax2.set_yticks([0.52, 0.54, 0.56, 0.58])
theme(ax2); sparse_year_ticks(ax2, top_row=True); tag(ax2, "B")

ax3.plot(X, ndep, color=PRIMARY, lw=1.55)
ax3.scatter(X, ndep, s=12, color=PRIMARY, zorder=3)
for xi, value in zip(X, ndep):
    ax3.text(xi, value + 0.65, f"{value:d}", fontsize=6.0,
             color=MUTED, ha="center", va="bottom")
ax3.set(xlim=(-0.25, 11.1), ylim=(40, 62), ylabel="Funded departments (n)")
theme(ax3); sparse_year_ticks(ax3); tag(ax3, "C")

rx = np.arange(len(ret))
ax4.plot(rx, ret, color=THIRD, lw=1.55)
ax4.scatter(rx, ret, s=12, color=THIRD, zorder=3)
ax4.text(9.18, ret[-1], f"{ret[-1]:.1f}%", color=THIRD,
         fontsize=6.9, va="center")
ax4.set(xlim=(-0.25, 11.1), ylim=(88, 100), ylabel="Retained in next FY (%)")
ax4.set_yticks([88, 92, 96, 100])
theme(ax4); sparse_year_ticks(ax4); tag(ax4, "D")
fig.subplots_adjust(left=0.10, right=0.955, bottom=0.10, top=0.97,
                    wspace=0.37, hspace=0.32)
save(fig, "figure2_concentration")


# FIGURE 3 -- paired classifier estimates by disease domain
s_a = A["subspecialty_real_funding"]
s_b = A["subspecialty_real_funding_B"]
order = ["Neuro-oncology", "Spine", "Cerebrovascular", "Functional/epilepsy",
         "Trauma/TBI", "Hydrocephalus", "Peripheral nerve"]
v_a = np.array([s_a[name] / 1e6 for name in order])
v_b = np.array([s_b[name] / 1e6 for name in order])
y = np.arange(len(order))
fig, ax = plt.subplots(figsize=(5.8, 3.7))
for yi, primary, robust in zip(y, v_a, v_b):
    ax.plot([robust, primary], [yi, yi], color=RULE, lw=2.0, zorder=1)
ax.scatter(v_b, y, s=27, facecolor="white", edgecolor=PRIMARY_L,
           linewidth=1.2, zorder=3)
ax.scatter(v_a, y, s=30, color=PRIMARY, zorder=4)
for yi, primary, robust in zip(y, v_a, v_b):
    ax.text(primary + 10, yi, f"{primary:.0f}", fontsize=6.8,
            color=PRIMARY, va="center")
    ax.text(robust - 10, yi, f"{robust:.0f}", fontsize=6.3,
            color=MUTED, ha="right", va="center")
ax.set_yticks(y); ax.set_yticklabels(order); ax.invert_yaxis()
ax.set(xlim=(0, 760), xlabel="Cumulative real funding, FY2015–FY2025 (US$ millions)")
theme(ax, grid_axis="x")
key = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=PRIMARY,
           markeredgecolor=PRIMARY, markersize=5, label="Expanded-text classifier"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
           markeredgecolor=PRIMARY_L, markersize=5,
           label="Title + spending-category classifier"),
]
ax.legend(handles=key, loc="lower right", bbox_to_anchor=(1.0, 1.005), ncol=2,
          frameon=False, handletextpad=0.4, columnspacing=1.2, borderaxespad=0)
fig.subplots_adjust(left=0.27, right=0.98, bottom=0.16, top=0.88)
save(fig, "figure3_domain")


# FIGURE 4 -- early-career transitions and T32 model estimates
pp = A["pi_pipeline"]
tm = pp["t32_adjusted_model"]
fig = plt.figure(figsize=(7.35, 3.15))
gs = fig.add_gridspec(1, 2, width_ratios=[0.9, 1.35], wspace=0.50)
ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[0, 1])
groups = ["Mentored K", "Fellowship\n(F30/F31/F32)"]
rates = np.array([pp["mentored_K_rate"], pp["fellowship_rate"]])
counts = [(pp["mentored_K_transitioned"], pp["mentored_K_cohort"]),
          (pp["fellowship_transitioned"], pp["fellowship_cohort"])]
yy = np.array([1, 0])
for yi, rate, color in zip(yy, rates, [PRIMARY, PRIMARY_L]):
    ax1.hlines(yi, 0, rate, color=RULE, lw=1.6, zorder=1)
    ax1.scatter(rate, yi, s=42, color=color, zorder=3)
for yi, rate, (num, den), color in zip(yy, rates, counts, [PRIMARY, PRIMARY_L]):
    ax1.text(rate + 1.4, yi, f"{rate:.0f}%  ({num}/{den})", color=color,
             fontsize=7.2, fontweight="bold", va="center")
ax1.set_yticks(yy); ax1.set_yticklabels(groups)
ax1.set(xlim=(0, 53), ylim=(-0.55, 1.55),
        xlabel="Later held an independent R-type award (%)")
ax1.set_xticks([0, 10, 20, 30, 40, 50])
theme(ax1, grid_axis="x"); tag(ax1, "A")

labels = ["Unadjusted", "Adjusted"]
irrs = np.array([tm["unadjusted_irr"], tm["adjusted_irr"]])
los = np.array([tm["unadjusted_irr_ci"][0], tm["adjusted_irr_ci"][0]])
his = np.array([tm["unadjusted_irr_ci"][1], tm["adjusted_irr_ci"][1]])
yy = np.array([1, 0])
colors = [FOURTH, PRIMARY]
ax2.axvline(1.0, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
for yi, irr, lo, hi, color in zip(yy, irrs, los, his, colors):
    ax2.plot([lo, hi], [yi, yi], color=color, lw=1.7, zorder=3)
    ax2.scatter(irr, yi, s=46, color=color, zorder=4)
    ax2.text(irr, yi + 0.17, f"{irr:.2f}  [{lo:.2f}–{hi:.2f}]", ha="center",
             va="bottom", color=INK, fontsize=6.9)
ax2.text(0.98, 1.035, "IRR  [95% CI]", transform=ax2.transAxes,
         ha="right", va="bottom", fontsize=6.7, color=MUTED)
ax2.set_xscale("log"); ax2.set(xlim=(0.28, 13.2), ylim=(-0.55, 1.55),
                              xlabel="Incidence rate ratio (distinct R01 PIs)")
major_ticks = [0.5, 1, 2, 4, 8]
ax2.xaxis.set_major_locator(FixedLocator(major_ticks))
ax2.xaxis.set_major_formatter(FixedFormatter(["0.5", "1", "2", "4", "8"]))
ax2.xaxis.set_minor_locator(NullLocator())
ax2.set_yticks(yy); ax2.set_yticklabels(labels)
theme(ax2, grid_axis="x"); tag(ax2, "B")
fig.subplots_adjust(left=0.12, right=0.99, bottom=0.19, top=0.93)
save(fig, "figure4_training")


# FIGURE 5 -- MeSH retrieval coverage
yrs = list(range(2015, 2026))
gen = np.array([M["years"][str(y)]["mesh_to_keyword_ratio"] for y in yrs])
jrn = np.array([M["years"][str(y)]["journal_mesh_ratio"] for y in yrs])
xi = np.arange(len(yrs))
fig, ax = plt.subplots(figsize=(6.15, 3.25))
ax.axvspan(8.5, 10.5, color=PALE, alpha=0.14, linewidth=0, zorder=0)
for values, color in ((gen, PRIMARY), (jrn, SECOND)):
    ax.plot(xi[:9], values[:9], color=color, lw=1.9, zorder=3)
    ax.plot(xi[8:], values[8:], color=color, lw=1.2, ls=(0, (3, 2)),
            alpha=0.75, zorder=3)
    ax.scatter(xi[:9], values[:9], s=11, color=color, zorder=4)
    ax.scatter(xi[9:], values[9:], s=12, facecolor="white", edgecolor=color,
               linewidth=0.9, zorder=4)
ax.text(9.5, 0.685, "incomplete\nindexing", ha="center", va="top",
        fontsize=6.6, color=MUTED, style="italic")
ax.text(10.18, gen[-1], f"General PubMed  {gen[-1]:.2f}", color=PRIMARY,
        fontsize=6.8, va="center")
ax.text(10.18, jrn[-1], f"Nine journals  {jrn[-1]:.2f}", color=SECOND,
        fontsize=6.8, va="center")
ax.set(xlim=(-0.35, 12.35), ylim=(0.20, 0.72), xlabel="Publication year",
       ylabel="MeSH retrieval coverage ratio")
ax.set_xticks(xi); ax.set_xticklabels([f"'{str(y)[2:]}" for y in yrs])
theme(ax)
fig.subplots_adjust(left=0.10, right=0.98, bottom=0.17, top=0.97)
save(fig, "figure5_mesh")


# SUPPLEMENTARY FIGURE 1 -- Lorenz curve, FY2025
if "fy2025_lorenz_curve" in A:
    xx = np.array(A["fy2025_lorenz_curve"]["department_share"], dtype=float)
    cum = np.array(A["fy2025_lorenz_curve"]["funding_share"], dtype=float)
else:
    cohort = pd.read_csv(os.path.join(OUT, "deduplicated_cohort.csv"))
    fy25 = cohort[cohort["Fiscal Year"] == 2025]
    inst = np.sort(fy25.groupby("Organization Name")["Total Cost"].sum().dropna().values)
    inst = inst[inst > 0]
    cum = np.insert(np.cumsum(inst) / inst.sum(), 0, 0)
    xx = np.linspace(0, 1, len(cum))
fig, ax = plt.subplots(figsize=(4.0, 3.85))
ax.plot([0, 1], [0, 1], color=PALE, lw=0.9, ls=(0, (3, 2)), zorder=1)
ax.plot(xx, cum, color=PRIMARY, lw=2.0, zorder=3)
ax.fill_between(xx, cum, xx, color=PRIMARY, alpha=0.08, linewidth=0, zorder=2)
ax.text(0.08, 0.86, f"Gini  {con['2025']['gini']:.3f}", fontsize=8.2,
        fontweight="bold", color=INK)
ax.text(0.60, 0.63, "line of equality", fontsize=6.4, color=MUTED,
        rotation=42, ha="center")
curve_y = np.interp(0.76, xx, cum)
ax.text(0.74, curve_y - 0.055, "FY2025 funding", fontsize=6.6,
        color=PRIMARY, rotation=52, ha="center")
ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Cumulative share of funded departments",
       ylabel="Cumulative share of funding")
ax.set_aspect("equal"); theme(ax, grid_axis="both")
fig.subplots_adjust(left=0.16, right=0.98, bottom=0.14, top=0.98)
save(fig, "sfigure1_lorenz")


# SUPPLEMENTARY FIGURE 2 -- NIH-wide benchmark, indexed real funding
dept = np.array([A["fy_real"][str(y)] for y in YEARS], dtype=float)
rpg = np.array([A["nih_rpg_real"][str(y)] for y in YEARS], dtype=float)
app = np.array([A["nih_appropriation_real"][str(y)] for y in YEARS], dtype=float)
dept_i, rpg_i, app_i = (100 * dept / dept[0], 100 * rpg / rpg[0],
                        100 * app / app[0])
fig, ax = plt.subplots(figsize=(6.2, 3.25))
series = [(dept_i, PRIMARY, "Neurosurgery cohort"),
          (rpg_i, SECOND, "NIH research grants"),
          (app_i, THIRD, "NIH appropriation")]
for values, color, label in series:
    ax.plot(X, values, color=color, lw=1.65, zorder=3)
    ax.scatter(X, values, s=10, color=color, zorder=4)
    ax.text(10.18, values[-1], f"{label}  {values[-1]:.0f}", color=color,
            fontsize=6.8, va="center")
ax.axhline(100, color=PALE, lw=0.8, ls=(0, (3, 2)))
ax.set(xlim=(-0.35, 12.8), ylim=(85, 210),
       ylabel="Real funding index (FY2015 = 100)")
sparse_year_ticks(ax); theme(ax)
fig.subplots_adjust(left=0.10, right=0.98, bottom=0.16, top=0.97)
save(fig, "sfigure2_benchmark")


# SUPPLEMENTARY FIGURE 3 -- deduplication sensitivity and yearly effect
ds = A["dedup_sensitivity"]
raw_real = np.array([ds["raw_real"][str(y)] / 1e6 for y in YEARS])
ddup_real = np.array([A["fy_real"][str(y)] / 1e6 for y in YEARS])
removed_pct = 100 * (raw_real - ddup_real) / raw_real
fig = plt.figure(figsize=(6.2, 4.05))
gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.0], hspace=0.08)
ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
ax1.fill_between(X, ddup_real, raw_real, color=PALE, alpha=0.15,
                 linewidth=0, zorder=1)
ax1.plot(X, raw_real, color=PALE, lw=1.45, zorder=3)
ax1.plot(X, ddup_real, color=PRIMARY, lw=1.85, zorder=4)
ax1.scatter(X, raw_real, s=9, color=PALE, zorder=5)
ax1.scatter(X, ddup_real, s=10, color=PRIMARY, zorder=5)
ax1.text(10.18, raw_real[-1], f"Before  ${raw_real[-1]:.0f}M", color=MUTED,
         fontsize=6.8, va="center")
ax1.text(10.18, ddup_real[-1], f"After  ${ddup_real[-1]:.0f}M", color=PRIMARY,
         fontsize=6.8, va="center")
ax1.set(xlim=(-0.35, 12.45), ylim=(120, 290),
        ylabel="Annual funding\n(2025 US$ millions)")
ax1.tick_params(axis="x", labelbottom=False); theme(ax1)
ax2.plot(X, removed_pct, color=INK, lw=1.2, zorder=3)
ax2.scatter(X, removed_pct, s=10, color=INK, zorder=4)
for xi, value in zip(X, removed_pct):
    ax2.text(xi, value + 0.65, f"{value:.1f}", fontsize=5.6,
             color=MUTED, ha="center", va="bottom")
ax2.set(ylim=(0, 15), ylabel="Removed\n(%)")
ax2.set_yticks([0, 5, 10, 15]); sparse_year_ticks(ax2); theme(ax2)
fig.subplots_adjust(left=0.11, right=0.98, bottom=0.13, top=0.98)
save(fig, "sfigure3_dedup")

print("\nAll figures written to", FIGDIR)
