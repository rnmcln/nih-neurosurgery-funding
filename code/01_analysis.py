#!/usr/bin/env python3
"""
01_analysis.py
==============
Reproducible analysis of NIH funding to neurosurgery departments,
FY1992-FY2025.

This script rebuilds every quantitative result from the raw NIH RePORTER
export. Key analytic safeguards include:
  * Records are DEDUPLICATED on Application ID (the prior pipeline summed
    192 duplicate export rows, modestly inflating funding totals).
  * CPI-U deflators are the official BLS annual averages (CUUR0000SA0),
    retrieved 2026-05-20; the 2025 annual average (321.943) is the BLS
    figure computed from 11 monthly observations (October 2025 was not
    collected because of the federal shutdown).
  * Adds Tier A hypothesis-testing analyses: PI-level career pipeline,
    institutional sustainability/survival, and formal trend statistics.

Input : a user-downloaded NIH RePORTER XLSX export (set REPORTER_FILE)
Output: outputs/analysis_output.json and outputs/deduplicated_cohort.csv
        (the outputs directory is intentionally excluded from version control)
"""

import json
import os
import re
import math
from collections import defaultdict, Counter

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# PATHS
# ----------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORTER_FILE = os.environ.get(
    "REPORTER_FILE",
    os.path.join(ROOT, "data", "source",
                 "SearchResult_Export_08Mar2026_022554.xlsx"),
)
OUT_DIR = os.environ.get("REPORTER_OUTPUT_DIR", os.path.join(ROOT, "outputs"))
os.makedirs(OUT_DIR, exist_ok=True)
if not os.path.exists(REPORTER_FILE):
    raise FileNotFoundError(
        "NIH RePORTER source export not found. Download the export described in "
        "DATA_PROVENANCE.md and set REPORTER_FILE to its local path."
    )

# ----------------------------------------------------------------------
# EXTERNAL REFERENCE CONSTANTS
# ----------------------------------------------------------------------
# CPI-U, US city average, all items, annual averages (BLS series CUUR0000SA0).
# Verified against the BLS public API on 2026-05-20.
CPI_U = {
    1992: 140.3, 1993: 144.5, 1994: 148.2, 1995: 152.4, 1996: 156.9,
    1997: 160.5, 1998: 163.0, 1999: 166.6, 2000: 172.2, 2001: 177.1,
    2002: 179.9, 2003: 184.0, 2004: 188.9, 2005: 195.3, 2006: 201.6,
    2007: 207.342, 2008: 215.303, 2009: 214.537, 2010: 218.056,
    2011: 224.939, 2012: 229.594, 2013: 232.957, 2014: 236.736,
    2015: 237.017, 2016: 240.007, 2017: 245.120, 2018: 251.107,
    2019: 255.657, 2020: 258.811, 2021: 270.970, 2022: 292.655,
    2023: 304.702, 2024: 313.689, 2025: 321.943,
}
CPI_BASE_YEAR = 2025  # constant 2025 dollars

# Total NIH appropriation (program level), USD millions, by fiscal year.
# Source: NIH Office of Budget / Congressional Research Service R43341.
# Used only as an approximate contextual denominator (see manuscript limitations).
NIH_APPROPRIATION_M = {
    2015: 30311, 2016: 32311, 2017: 34301, 2018: 37311, 2019: 39313,
    2020: 41690, 2021: 42937, 2022: 45180, 2023: 47684, 2024: 47099,
    2025: 47099,
}

# NIH-wide Research Project Grant (RPG) obligations and award counts.
# Source: NIH Data Book, Research Project Grants, FY2015-FY2025
# (https://report.nih.gov/nihdatabook/category/4; accessed 2026-10-12).
NIH_RPG_FUNDING_USD = {
    2015: 15862012059, 2016: 17137754907, 2017: 18321187243,
    2018: 19869684598, 2019: 21589062208, 2020: 22636319154,
    2021: 23280221639, 2022: 24399835832, 2023: 25365442608,
    2024: 25721050507, 2025: 26503031656,
}
NIH_RPG_AWARD_COUNTS = {2023: 41597, 2024: 41470, 2025: 39239}

ANALYSIS_START, ANALYSIS_END = 2015, 2025  # primary analysis window

# ----------------------------------------------------------------------
# LOAD + DEDUPLICATE
# ----------------------------------------------------------------------
print("=" * 72)
print("LOADING NIH RePORTER EXPORT")
print("=" * 72)
df = pd.read_excel(REPORTER_FILE)
n_raw = len(df)
print(f"Raw export rows                : {n_raw}")

# Total Cost: coerce, treat missing as 0 for summation transparency.
df["Total Cost"] = pd.to_numeric(df["Total Cost"], errors="coerce")
df["cost"] = df["Total Cost"].fillna(0.0)
df["Fiscal Year"] = pd.to_numeric(df["Fiscal Year"], errors="coerce").astype("Int64")
df["Activity"] = df["Activity"].fillna("").astype(str).str.strip().str.upper()

# Deduplicate on Application ID. RePORTER exports can repeat application
# records; the prior pipeline summed all rows. We keep, per Application ID,
# the row with the largest Total Cost (resolves exact duplicates and the
# null-vs-value case identically).
dup_app = df["Application ID"].duplicated().sum()
df_sorted = df.sort_values("cost", ascending=False)
dedup = df_sorted.drop_duplicates(subset=["Application ID"], keep="first").copy()
n_dedup = len(dedup)
print(f"Duplicate Application ID rows  : {dup_app}")
print(f"Deduplicated grant-year records: {n_dedup}")
print(f"Unique projects (Project Number): {dedup['Project Number'].nunique()}")
print(f"Unique institutions            : {dedup['Organization Name'].nunique()}")
print(f"Unique contact PIs             : {dedup['Contact PI Person ID'].nunique()}")

raw_total = df["cost"].sum()
dedup_total = dedup["cost"].sum()
print(f"Total funding, all rows        : ${raw_total:,.0f}")
print(f"Total funding, deduplicated    : ${dedup_total:,.0f}")
print(f"Double-counting removed        : ${raw_total - dedup_total:,.0f} "
      f"({100*(raw_total-dedup_total)/raw_total:.2f}%)")

# ----------------------------------------------------------------------
# SUBSPECIALTY CLASSIFIER  (defined up front; applied to `dedup`)
# ----------------------------------------------------------------------
# Multi-label keyword screen. Long/multi-word terms use substring matching;
# short abbreviations use word-boundary matching to avoid false positives.
SUBSPECIALTY_KEYWORDS = {
    "Neuro-oncology": {
        "phrases": ["glioma", "glioblastoma", "brain tumor", "brain neoplasm",
                    "meningioma", "brain cancer", "neuro-oncol", "neurooncol",
                    "medulloblastoma", "craniopharyngioma", "schwannoma",
                    "pituitary", "brain metastas", "ependymoma", "astrocytoma"],
        "abbrev": ["gbm"],
    },
    "Cerebrovascular": {
        "phrases": ["aneurysm", "cerebrovascular", "stroke",
                    "subarachnoid hemorrhage", "arteriovenous malformation",
                    "cavernous malformation", "intracranial hemorrhage",
                    "cerebral hemorrhage", "moyamoya", "vasospasm", "endovascular",
                    "carotid"],
        "abbrev": ["sah", "avm", "ich"],
    },
    "Spine": {
        "phrases": ["spinal", "spine", "lumbar", "cervical", "disc herniation",
                    "scoliosis", "vertebr", "spinal cord injury", "spinal cord",
                    "laminectomy", "degenerative disc", "myelopathy",
                    "radiculopathy", "intervertebral"],
        "abbrev": ["sci"],
    },
    "Trauma/TBI": {
        "phrases": ["traumatic brain injury", "craniocerebral trauma", "subdural",
                    "epidural hematoma", "skull fracture", "concussion",
                    "head injury", "blast injury", "neurotrauma",
                    "decompressive craniectomy"],
        "abbrev": ["tbi"],
    },
    "Functional/epilepsy": {
        "phrases": ["deep brain stimulation", "epilepsy", "seizure",
                    "functional neurosurg", "movement disorder", "parkinson",
                    "tremor", "dystonia", "responsive neurostimulation",
                    "vagus nerve stimulat", "stereoelectroencephalography"],
        "abbrev": ["dbs"],
    },
    "Hydrocephalus": {
        "phrases": ["hydrocephalus", "ventriculoperitoneal", "cerebrospinal fluid",
                    "normal pressure hydrocephalus",
                    "endoscopic third ventriculostomy"],
        "abbrev": ["nph", "csf"],
    },
    "Peripheral nerve": {
        "phrases": ["peripheral nerve", "brachial plexus", "nerve injury",
                    "nerve repair", "nerve regeneration", "nerve graft",
                    "nerve transfer", "carpal tunnel", "ulnar nerve",
                    "median nerve", "radial nerve", "sciatic nerve",
                    "neurofibromatosis", "neurofibroma",
                    "peripheral neuropath"],
        "abbrev": [],
    },
}
STRATEGY_FIELDS = {
    "A_fulltext": ["Project Title", "Project Terms", "Project Abstract",
                   "NIH Spending Categorization", "Public Health Relevance"],
    "B_title_spend": ["Project Title", "NIH Spending Categorization"],
}


def _match(text, spec):
    for p in spec["phrases"]:
        if p in text:
            return True
    for a in spec["abbrev"]:
        if re.search(r"\b" + re.escape(a) + r"\b", text):
            return True
    return False


def classify(row, fields):
    text = " ".join(str(row.get(f, "") or "") for f in fields).lower()
    labels = [name for name, spec in SUBSPECIALTY_KEYWORDS.items()
              if _match(text, spec)]
    return set(labels) if labels else {"Other/Unclassified"}


dedup["subspecialties"] = dedup.apply(
    lambda r: classify(r, STRATEGY_FIELDS["A_fulltext"]), axis=1)
dedup["subspecialties_B"] = dedup.apply(
    lambda r: classify(r, STRATEGY_FIELDS["B_title_spend"]), axis=1)

# ----------------------------------------------------------------------
# 1. FUNDING BY FISCAL YEAR  (nominal + constant 2025 USD)
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("1. FUNDING BY FISCAL YEAR")
print("=" * 72)
fy_nominal = dedup.groupby("Fiscal Year")["cost"].sum().to_dict()
fy_nominal = {int(k): float(v) for k, v in fy_nominal.items() if pd.notna(k)}
fy_records = dedup.groupby("Fiscal Year").size().to_dict()
fy_records = {int(k): int(v) for k, v in fy_records.items() if pd.notna(k)}
fy_projects = dedup.groupby("Fiscal Year")["Project Number"].nunique().to_dict()
fy_projects = {int(k): int(v) for k, v in fy_projects.items() if pd.notna(k)}


def to_real(nominal, year):
    """Deflate a nominal USD amount in `year` to constant CPI_BASE_YEAR dollars."""
    return nominal * CPI_U[CPI_BASE_YEAR] / CPI_U[year]


fy_real = {y: to_real(v, y) for y, v in fy_nominal.items() if y in CPI_U}

for y in range(1992, 2027):
    if y in fy_nominal:
        r = fy_real.get(y, float("nan"))
        print(f"  FY{y}: nominal ${fy_nominal[y]:>13,.0f}   real ${r:>13,.0f}   "
              f"records {fy_records.get(y,0):>4}   projects {fy_projects.get(y,0):>4}")


def cagr(v0, v1, n):
    return (v1 / v0) ** (1.0 / n) - 1.0


nom_cagr = cagr(fy_nominal[2015], fy_nominal[2025], 10)
real_cagr = cagr(fy_real[2015], fy_real[2025], 10)
print(f"\n  Nominal CAGR FY2015-FY2025: {100*nom_cagr:.2f}%")
print(f"  Real CAGR    FY2015-FY2025: {100*real_cagr:.2f}%")
print(f"  Real funding FY2015: ${fy_real[2015]:,.0f}")
print(f"  Real funding FY2023: ${fy_real[2023]:,.0f}  (peak)")
print(f"  Real funding FY2025: ${fy_real[2025]:,.0f}")
print(f"  Peak-to-FY2025 real change: {100*(fy_real[2025]/fy_real[2023]-1):.2f}%")

# ----------------------------------------------------------------------
# 1b. DEDUPLICATION SENSITIVITY
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("1b. DEDUPLICATION SENSITIVITY ANALYSIS")
print("=" * 72)
# Raw (pre-deduplication) annual funding, all 8,672 export rows.
raw_fy_nominal = df.groupby("Fiscal Year")["cost"].sum().to_dict()
raw_fy_nominal = {int(k): float(v) for k, v in raw_fy_nominal.items() if pd.notna(k)}
raw_fy_real = {y: to_real(v, y) for y, v in raw_fy_nominal.items() if y in CPI_U}
raw_nom_cagr = cagr(raw_fy_nominal[2015], raw_fy_nominal[2025], 10)
raw_real_cagr = cagr(raw_fy_real[2015], raw_fy_real[2025], 10)
print(f"  {'FY':>6} {'raw nominal':>15} {'dedup nominal':>15} {'diff %':>8}")
dedup_sens = {}
for y in range(ANALYSIS_START, ANALYSIS_END + 1):
    rn, dn = raw_fy_nominal[y], fy_nominal[y]
    dedup_sens[y] = {"raw_nominal": rn, "dedup_nominal": dn,
                     "raw_real": raw_fy_real[y], "dedup_real": fy_real[y]}
    print(f"  FY{y}: {rn:>15,.0f} {dn:>15,.0f} {100*(rn-dn)/rn:>7.2f}%")
print(f"  Raw real CAGR FY2015-2025   : {100*raw_real_cagr:.2f}%  "
      f"(deduplicated: {100*real_cagr:.2f}%)")
# Does the FY2023 peak / FY2024-25 decline survive in the raw series?
raw_peak_year = max(range(2015, 2026), key=lambda y: raw_fy_real[y])
print(f"  Raw-series real-funding peak year: FY{raw_peak_year} "
      f"(deduplicated peak: FY{max(range(2015,2026), key=lambda y: fy_real[y])})")
print(f"  Raw FY2023->FY2025 real change: "
      f"{100*(raw_fy_real[2025]/raw_fy_real[2023]-1):.2f}%  "
      f"(deduplicated: {100*(fy_real[2025]/fy_real[2023]-1):.2f}%)")

# Characterise the duplicate Application-ID rows.
dup_mask = df["Application ID"].duplicated(keep=False)
dup_rows = df[dup_mask].copy()
n_dup_appids = dup_rows["Application ID"].nunique()
exact_dups = 0
costdiff_dups = 0
for aid, grp in dup_rows.groupby("Application ID"):
    costs = grp["Total Cost"]
    nonnull = costs.dropna().unique()
    if costs.isna().any() and len(nonnull) <= 1:
        costdiff_dups += 1            # one copy missing a cost value
    elif len(set(costs.fillna(-1).round(2))) == 1:
        exact_dups += 1               # identical cost across copies
    else:
        costdiff_dups += 1
dup_type_dist = dup_rows["Type"].value_counts(dropna=False).to_dict()
print(f"\n  Duplicate Application IDs: {n_dup_appids} "
      f"(involving {len(dup_rows)} export rows)")
print(f"   - copies with identical Total Cost   : {exact_dups}")
print(f"   - copies differing only by a missing/updated cost value: {costdiff_dups}")
print(f"   Type-code distribution among duplicate rows: "
      f"{ {str(k): int(v) for k, v in dup_type_dist.items()} }")

# ----------------------------------------------------------------------
# 2. BUDGET SHARE  (nominal departmental funding / total NIH appropriation)
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("2. SHARE OF TOTAL NIH APPROPRIATION  (approximate benchmark)")
print("=" * 72)
budget_share = {}
for y in range(ANALYSIS_START, ANALYSIS_END + 1):
    share = 100.0 * fy_nominal[y] / (NIH_APPROPRIATION_M[y] * 1e6)
    budget_share[y] = share
    print(f"  FY{y}: ${fy_nominal[y]/1e6:7.1f}M / ${NIH_APPROPRIATION_M[y]/1000:.1f}B "
          f"= {share:.3f}%")

# Contextualize the post-FY2023 departmental decline against two NIH-wide
# series with different denominators: total NIH appropriations and RPG
# obligations. Comparisons are indexed within series and inflation-adjusted.
nih_rpg_real = {y: to_real(v, y) for y, v in NIH_RPG_FUNDING_USD.items()}
nih_appropriation_real = {
    y: to_real(v * 1e6, y) for y, v in NIH_APPROPRIATION_M.items()
}
nih_benchmark = {
    "department_real_change_2023_2025_pct":
        float(100 * (fy_real[2025] / fy_real[2023] - 1)),
    "nih_rpg_real_change_2023_2025_pct":
        float(100 * (nih_rpg_real[2025] / nih_rpg_real[2023] - 1)),
    "nih_appropriation_real_change_2023_2025_pct":
        float(100 * (nih_appropriation_real[2025] /
                     nih_appropriation_real[2023] - 1)),
    "nih_rpg_award_count_change_2023_2025_pct":
        float(100 * (NIH_RPG_AWARD_COUNTS[2025] /
                     NIH_RPG_AWARD_COUNTS[2023] - 1)),
}
print("\n  Post-FY2023 real-dollar benchmark (FY2023 to FY2025):")
print(f"    Neurosurgery-department cohort: "
      f"{nih_benchmark['department_real_change_2023_2025_pct']:.1f}%")
print(f"    NIH-wide RPG obligations:       "
      f"{nih_benchmark['nih_rpg_real_change_2023_2025_pct']:.1f}%")
print(f"    Total NIH appropriation:        "
      f"{nih_benchmark['nih_appropriation_real_change_2023_2025_pct']:.1f}%")
print(f"    NIH-wide RPG award count:       "
      f"{nih_benchmark['nih_rpg_award_count_change_2023_2025_pct']:.1f}%")

# ----------------------------------------------------------------------
# 3. INSTITUTIONAL CONCENTRATION  (Gini, HHI, top-k shares)
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("3. INSTITUTIONAL CONCENTRATION, FY2015-FY2025")
print("=" * 72)


def gini(values):
    """Gini coefficient of a list of non-negative values (0=equal, 1=maximal)."""
    x = np.sort(np.asarray([v for v in values if v >= 0], dtype=float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return float("nan")
    cum = np.cumsum(x)
    # mean absolute difference formulation
    return (2.0 * np.sum((np.arange(1, n + 1)) * x) / (n * cum[-1])) - (n + 1.0) / n


def hhi(values):
    """Herfindahl-Hirschman Index on shares (sum of squared percentage shares)."""
    x = np.asarray([v for v in values if v > 0], dtype=float)
    if x.sum() == 0:
        return float("nan")
    s = 100.0 * x / x.sum()
    return float(np.sum(s ** 2))


concentration = {}
for y in range(ANALYSIS_START, ANALYSIS_END + 1):
    sub = dedup[dedup["Fiscal Year"] == y]
    inst = sub.groupby("Organization Name")["cost"].sum()
    inst = inst[inst > 0].sort_values(ascending=False)
    tot = inst.sum()
    g = gini(inst.values)
    h = hhi(inst.values)
    top5 = 100.0 * inst.iloc[:5].sum() / tot
    top10 = 100.0 * inst.iloc[:10].sum() / tot
    cv = float(np.std(inst.values, ddof=1) / np.mean(inst.values))
    concentration[y] = {
        "n_departments": int(len(inst)),
        "gini": float(g),
        "hhi": float(h),
        "top5_pct": float(top5),
        "top10_pct": float(top10),
        "cv": cv,
        "median": float(inst.median()),
        "iqr_low": float(inst.quantile(0.25)),
        "iqr_high": float(inst.quantile(0.75)),
    }
    print(f"  FY{y}: n={len(inst):2d}  Gini={g:.3f}  HHI={h:6.0f}  "
          f"top5={top5:5.1f}%  top10={top10:5.1f}%  median=${inst.median():,.0f}")

# Rights-safe, anonymous aggregate coordinates for the FY2025 Lorenz curve.
# These shares allow the public figure package to reproduce Figure S1 without
# redistributing institution names, institution-level amounts, or row-level data.
fy2025_inst = (dedup[dedup["Fiscal Year"] == 2025]
               .groupby("Organization Name")["cost"].sum())
fy2025_inst = np.sort(fy2025_inst[fy2025_inst > 0].values)
fy2025_lorenz_funding = np.insert(
    np.cumsum(fy2025_inst) / fy2025_inst.sum(), 0, 0)
fy2025_lorenz_departments = np.linspace(0, 1, len(fy2025_lorenz_funding))

# Cumulative top institutions FY2015-2025 (real USD)
dedup_win = dedup[(dedup["Fiscal Year"] >= ANALYSIS_START) &
                  (dedup["Fiscal Year"] <= ANALYSIS_END)].copy()
dedup_win["real_cost"] = dedup_win.apply(
    lambda r: to_real(r["cost"], int(r["Fiscal Year"])), axis=1)
inst_cum = dedup_win.groupby("Organization Name")["real_cost"].sum().sort_values(ascending=False)
print("\n  Top 10 institutions, cumulative real funding FY2015-FY2025:")
for i, (org, val) in enumerate(inst_cum.head(10).items(), 1):
    print(f"   {i:2d}. {org[:45]:45s} ${val:,.0f}")

# ----------------------------------------------------------------------
# 4. MECHANISM MIX  (activity codes)
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("4. FUNDING MECHANISM MIX, FY2015-FY2025 (real USD)")
print("=" * 72)
mech = dedup_win.groupby("Activity")["real_cost"].agg(["sum", "count"]).sort_values(
    "sum", ascending=False)
win_total_real = dedup_win["real_cost"].sum()
mech_mix = {}
for act, row in mech.head(15).iterrows():
    pct = 100.0 * row["sum"] / win_total_real
    mech_mix[act] = {"real_funding": float(row["sum"]), "n": int(row["count"]),
                     "pct": float(pct)}
    print(f"  {act:6s}: ${row['sum']:>13,.0f}  ({pct:5.1f}%)  n={int(row['count'])}")

r01_r37 = mech.loc[mech.index.isin(["R01", "R37"]), "sum"].sum()
r01_only = mech.loc["R01", "sum"] if "R01" in mech.index else 0
print(f"\n  R01+R37 share: {100*r01_r37/win_total_real:.1f}%")
print(f"  R01 only share: {100*r01_only/win_total_real:.1f}%")

# Center grants P01/P30/P50 share, first vs last year
center = ["P01", "P30", "P50"]
for y in (2015, 2025):
    sub = dedup[dedup["Fiscal Year"] == y]
    c = sub[sub["Activity"].isin(center)]["cost"].sum()
    t = sub["cost"].sum()
    print(f"  Center grants (P01/P30/P50) FY{y}: {100*c/t:.1f}% of funding")

# ----------------------------------------------------------------------
# 5. SUBSPECIALTY FUNDING  (multi-label keyword screen + robustness check)
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("5. SUBSPECIALTY FUNDING, FY2015-FY2025")
print("=" * 72)

# --- Subspecialty funding under strategy A (full text), real USD ---
sub_real = defaultdict(float)
sub_count = defaultdict(int)
sub_fy_real = defaultdict(lambda: defaultdict(float))
for _, row in dedup_win.iterrows():
    for s in row["subspecialties"]:
        sub_real[s] += row["real_cost"]
        sub_count[s] += 1
        sub_fy_real[s][int(row["Fiscal Year"])] += row["real_cost"]
print("  Cumulative real funding by subspecialty (strategy A, full text):")
for s, v in sorted(sub_real.items(), key=lambda x: -x[1]):
    print(f"    {s:22s} ${v:>13,.0f}  (n={sub_count[s]})")

# --- Robustness: independent strategy B (title + spending category only) ---
sub_real_B = defaultdict(float)
for _, row in dedup_win.iterrows():
    for s in row["subspecialties_B"]:
        sub_real_B[s] += row["real_cost"]
print("\n  Robustness check -- cumulative real funding (strategy B, "
      "title+spending only):")
order_A = [s for s, _ in sorted(sub_real.items(), key=lambda x: -x[1])
           if s != "Other/Unclassified"]
order_B = [s for s, _ in sorted(sub_real_B.items(), key=lambda x: -x[1])
           if s != "Other/Unclassified"]
for s in order_B:
    print(f"    {s:22s} ${sub_real_B[s]:>13,.0f}")

# Concordance metrics between the two independent classifications
try:
    from scipy.stats import spearmanr as _spear
    cats6 = [c for c in SUBSPECIALTY_KEYWORDS]
    rankA = [order_A.index(c) for c in cats6]
    rankB = [order_B.index(c) for c in cats6]
    rho_rank, p_rank = _spear(rankA, rankB)
except Exception:
    rho_rank = p_rank = None

# Per-grant agreement (Jaccard of label sets, A vs B) over the window
jacc = []
exact_match = 0
for _, row in dedup_win.iterrows():
    a, b = row["subspecialties"], row["subspecialties_B"]
    inter = len(a & b)
    union = len(a | b)
    jacc.append(inter / union if union else 1.0)
    exact_match += (a == b)
mean_jaccard = float(np.mean(jacc))
exact_pct = 100.0 * exact_match / len(dedup_win)
unclassified_A = 100.0 * sum(
    dedup_win["subspecialties"].apply(lambda s: s == {"Other/Unclassified"})) / len(dedup_win)
print(f"\n  Strategy A vs B: rank-order Spearman rho = "
      f"{rho_rank if rho_rank is None else round(rho_rank,3)}")
print(f"  Mean per-grant Jaccard (A vs B label sets) = {mean_jaccard:.3f}")
print(f"  Exact label-set agreement (A vs B) = {exact_pct:.1f}%")
print(f"  Grants unclassified under strategy A = {unclassified_A:.1f}%")

# ----------------------------------------------------------------------
# 6b. CROSS-WALK: count of grants per subspecialty (for transparency)
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# 6. TIER A  --  PI-LEVEL CAREER PIPELINE
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("6. TIER A: SURGEON-SCIENTIST / PI CAREER PIPELINE")
print("=" * 72)

FELLOWSHIP = {"F30", "F31", "F32", "F99"}
MENTORED_K = {"K01", "K07", "K08", "K12", "K22", "K23", "K25", "K38", "K99", "K00"}
EARLY_CAREER = FELLOWSHIP | MENTORED_K | {"T32"}
INDEPENDENT_R = {"R01", "R37", "RF1", "R35", "DP1", "DP2", "DP5", "R56", "R29"}

# Per PI: list of (FY, activity)
pi_hist = defaultdict(list)
for _, row in dedup.iterrows():
    pid = row["Contact PI Person ID"]
    fy = row["Fiscal Year"]
    act = row["Activity"]
    if pd.notna(pid) and pd.notna(fy):
        pi_hist[pid].append((int(fy), act))

pi_first_career = {}
pi_first_indep = {}
pi_first_any = {}
for pid, hist in pi_hist.items():
    hist_sorted = sorted(hist, key=lambda item: (item[0], str(item[1])))
    pi_first_any[pid] = hist_sorted[0][0]
    cyears = [fy for fy, a in hist_sorted if a in (FELLOWSHIP | MENTORED_K)]
    ryears = [fy for fy, a in hist_sorted if a in INDEPENDENT_R]
    if cyears:
        pi_first_career[pid] = min(cyears)
    if ryears:
        pi_first_indep[pid] = min(ryears)

# Cohort: PIs whose first NIH award in this dataset is a mentored-K or
# fellowship award (entered as early-career within a neurosurgery department).
# Headline cohort requires first award <= FY2019 (>=6 years potential
# follow-up, exceeding the observed median time-to-R01); a wider cohort
# (<= FY2022) is reported with an explicit right-censoring caveat.


def pipeline_stats(cutoff):
    cohort = [pid for pid, fy0 in pi_first_career.items()
              if sorted(pi_hist[pid], key=lambda item: (item[0], str(item[1])))[0][1]
              in (FELLOWSHIP | MENTORED_K)
              and fy0 <= cutoff]
    trans, tt = [], []
    for pid in cohort:
        fy0 = pi_first_career[pid]
        if pid in pi_first_indep and pi_first_indep[pid] >= fy0:
            trans.append(pid)
            tt.append(pi_first_indep[pid] - fy0)
    return cohort, trans, tt


pipeline_cohort, transitioned, times = pipeline_stats(2019)   # mature cohort
wide_cohort, wide_trans, _ = pipeline_stats(2022)             # wider, censored
pipeline_sensitivity = {}
for cutoff in (2017, 2018, 2019, 2020):
    c, t, tt = pipeline_stats(cutoff)
    pipeline_sensitivity[cutoff] = {
        "cohort_n": len(c),
        "transitioned_n": len(t),
        "transition_rate": float(100 * len(t) / len(c)) if c else None,
        "median_years_to_R": float(np.median(tt)) if tt else None,
    }

n_cohort = len(pipeline_cohort)
n_trans = len(transitioned)
print(f"  Mature early-career-entry PI cohort (first award = F/K, FY<=2019): "
      f"n={n_cohort}")
print(f"  Subsequently held an independent R award (R01/R37/RF1/...): "
      f"n={n_trans} ({100*n_trans/n_cohort:.1f}%)")
if times:
    print(f"  Median years to first independent R award: {np.median(times):.1f} "
          f"(IQR {np.percentile(times,25):.1f}-{np.percentile(times,75):.1f})")
print(f"  Wider cohort (first award FY<=2022, right-censored): n={len(wide_cohort)}, "
      f"transitioned n={len(wide_trans)} ({100*len(wide_trans)/len(wide_cohort):.1f}%)")
print("  Cutoff sensitivity (transitioned/cohort): " + ", ".join(
    f"FY<={y}: {v['transitioned_n']}/{v['cohort_n']} "
    f"({v['transition_rate']:.1f}%)"
    for y, v in pipeline_sensitivity.items()))

# Mentored-K subgroup specifically (K08/K23 are the classic clinician path)
k_cohort = [pid for pid in pipeline_cohort
            if sorted(pi_hist[pid], key=lambda item: (item[0], str(item[1])))[0][1]
            in MENTORED_K]
k_trans = [pid for pid in k_cohort if pid in transitioned]
print(f"  Mentored-K entrants (mature cohort): n={len(k_cohort)}; "
      f"transitioned to R: n={len(k_trans)} "
      f"({100*len(k_trans)/max(1,len(k_cohort)):.1f}%)")
f_cohort = [pid for pid in pipeline_cohort
            if sorted(pi_hist[pid], key=lambda item: (item[0], str(item[1])))[0][1]
            in FELLOWSHIP]
f_trans = [pid for pid in f_cohort if pid in transitioned]
print(f"  Fellowship (F) entrants (mature cohort): n={len(f_cohort)}; "
      f"transitioned to R: n={len(f_trans)} "
      f"({100*len(f_trans)/max(1,len(f_cohort)):.1f}%)")

# Multi-PI structure. RePORTER provides a Person ID only for the contact PI;
# additional PIs appear as names in "Other PI or Project Leader(s)". The
# pipeline analysis therefore tracks contact-PI status only.
opl = dedup["Other PI or Project Leader(s)"].astype(str).str.strip()
has_multi_pi = ~opl.isin(["", "nan", "None", "Not Applicable"])
multi_pi_pct = 100.0 * has_multi_pi.sum() / len(dedup)
print(f"  Grant-years with >1 listed PI (multi-PI awards): "
      f"{has_multi_pi.sum()} ({multi_pi_pct:.1f}%); "
      f"pipeline analysis tracks contact PI only.")

# PI trajectories follow the stable Contact PI Person ID across institutions;
# institution-level concentration assigns each grant-year to the organization
# named on that record. Quantify observed institutional mobility explicitly.
pi_org_counts_all = (dedup.dropna(subset=["Contact PI Person ID",
                                          "Organization Name"])
                     .groupby("Contact PI Person ID")["Organization Name"]
                     .nunique())
pi_org_counts_win = (dedup_win.dropna(subset=["Contact PI Person ID",
                                              "Organization Name"])
                     .groupby("Contact PI Person ID")["Organization Name"]
                     .nunique())
pi_mobility = {
    "all_years_multi_institution_n": int((pi_org_counts_all > 1).sum()),
    "all_years_total_pi_n": int(len(pi_org_counts_all)),
    "fy2015_2025_multi_institution_n": int((pi_org_counts_win > 1).sum()),
    "fy2015_2025_total_pi_n": int(len(pi_org_counts_win)),
    "trajectory_attribution": "Contact PI Person ID across organizations",
    "institution_attribution": "Organization Name on each grant-year record",
}
print(f"  PIs linked to >1 institution in FY2015-FY2025: "
      f"{pi_mobility['fy2015_2025_multi_institution_n']}/"
      f"{pi_mobility['fy2015_2025_total_pi_n']}; PI trajectories follow "
      "Person ID, while institutional metrics use the organization on each record.")

# Focused manual audit of the 20 highest-funded contact PIs (nominal funding,
# FY2015-FY2025). Categories were verified against official institutional
# profiles on 2026-10-12; no biographical text is reproduced.
PI_TYPE_AUDIT = {
    3160894: ("neurosurgeon", "https://neurosurgery.duke.edu/john-h-sampson"),
    9729256: ("neurosurgeon", "https://neuroscience.ufl.edu/profile/hoh-brian/"),
    7351737: ("neurosurgeon", "https://www.feinberg.northwestern.edu/faculty-profiles/az/profile.html?xid=34049"),
    8926441: ("neurosurgeon", "https://braintumorcenter.ucsf.edu/people/edward-chang"),
    1896710: ("non-physician scientist", "https://med.stanford.edu/profiles/ivan-soltesz"),
    7690446: ("non-neurosurgeon physician", "https://www.ucsfhealth.org/providers/susan-chang"),
    1894343: ("non-neurosurgeon physician", "https://www.med.upenn.edu/apps/faculty/index.php/g20000440/p10366"),
    1898964: ("neurosurgeon", "https://braintumorcenter.ucsf.edu/people/mitchel-berger"),
    2190844: ("neurosurgeon", "https://faculty.mdanderson.org/profiles/frederick_lang.html"),
    9424298: ("neurosurgeon", "https://www.bcm.edu/people-search/sameer-sheth-30585"),
    1900873: ("non-neurosurgeon physician", "https://neuroscience.ufl.edu/profile/mitchell-duane/"),
    3143458: ("neurosurgeon", "https://www.uclahealth.org/providers/linda-liau"),
    1887100: ("non-physician scientist", "https://braintumorcenter.ucsf.edu/people/arturo-alvarez-buylla"),
    1989991: ("neurosurgeon", "https://medicine.osu.edu/find-a-researcher/krzysztof-bankiewicz-100004112"),
    8544696: ("neurosurgeon", "https://profiles.stanford.edu/jaimie-henderson"),
    1878386: ("neurosurgeon", "https://medicine.yale.edu/profile/steven-schiff/"),
    6278770: ("neurosurgeon", "https://profiles.utsouthwestern.edu/profile/205161/nader-pouratian.html"),
    6066228: ("neurosurgeon", "https://neurosurgery.wustl.edu/people/gregory-j-zipfel-md/"),
    1928282: ("neurosurgeon", "https://med.emory.edu/departments/neurosurgery/education/academic_tracks/functional_track/components.html"),
    1865771: ("neurosurgeon", "https://med.stanford.edu/profiles/gary-steinberg"),
}
top_pi = (dedup_win.dropna(subset=["Contact PI Person ID"])
          .groupby(["Contact PI Person ID", "Contact PI / Project Leader"],
                   dropna=False)["cost"].sum()
          .sort_values(ascending=False).head(20).reset_index())
top_pi_audit = []
for rank, row in top_pi.iterrows():
    pid = int(row["Contact PI Person ID"])
    category, source_url = PI_TYPE_AUDIT[pid]
    top_pi_audit.append({
        "rank": int(rank + 1),
        "contact_pi_person_id": pid,
        "contact_pi_name": str(row["Contact PI / Project Leader"]),
        "nominal_funding_fy2015_2025": float(row["cost"]),
        "category": category,
        "source_url": source_url,
    })
top_pi_category_counts = Counter(r["category"] for r in top_pi_audit)
print("  Manual audit, 20 highest-funded contact PIs: " + ", ".join(
    f"{category} n={count}"
    for category, count in sorted(top_pi_category_counts.items())))

# Institution-level: does T32 presence predict R01 portfolio size?
inst_t32_years = defaultdict(set)
inst_r01_pis_late = defaultdict(set)
inst_r01_funding_late = defaultdict(float)
for _, row in dedup.iterrows():
    org = row["Organization Name"]
    fy = row["Fiscal Year"]
    act = row["Activity"]
    if pd.isna(fy):
        continue
    fy = int(fy)
    if act == "T32" and 2010 <= fy <= 2020:
        inst_t32_years[org].add(fy)
    if act == "R01" and 2021 <= fy <= 2025:
        inst_r01_pis_late[org].add(row["Contact PI Person ID"])
        inst_r01_funding_late[org] += row["cost"]

all_orgs = set(dedup["Organization Name"].dropna())
t32_orgs = [o for o in all_orgs if o in inst_t32_years]
no_t32_orgs = [o for o in all_orgs if o not in inst_t32_years]
r01_t32 = [len(inst_r01_pis_late.get(o, set())) for o in t32_orgs]
r01_no = [len(inst_r01_pis_late.get(o, set())) for o in no_t32_orgs]
print(f"\n  Institution-level T32 -> R01 association:")
print(f"    Institutions with a T32 (FY2010-2020): n={len(t32_orgs)}; "
      f"median distinct R01 PIs FY2021-2025 = {np.median(r01_t32):.1f}")
print(f"    Institutions without a T32:            n={len(no_t32_orgs)}; "
      f"median distinct R01 PIs FY2021-2025 = {np.median(r01_no):.1f}")
try:
    from scipy.stats import mannwhitneyu, spearmanr
    u, p_u = mannwhitneyu(r01_t32, r01_no, alternative="greater")
    print(f"    Mann-Whitney U (T32 > no-T32), one-sided p = {p_u:.4g}")
    # dose-response: cumulative T32 years vs R01 PI count
    xs = [len(inst_t32_years.get(o, set())) for o in all_orgs]
    ys = [len(inst_r01_pis_late.get(o, set())) for o in all_orgs]
    rho, p_rho = spearmanr(xs, ys)
    print(f"    Spearman rho (T32 years vs R01 PI count) = {rho:.3f}, p = {p_rho:.4g}")
except Exception as e:  # pragma: no cover
    p_u = rho = p_rho = None
    print(f"    [scipy unavailable: {e}]")

# --- Adjusted model: does the T32 association survive control for baseline
#     departmental research intensity? ---
# Baseline (FY2010-FY2020) total NIH funding and distinct R01-PI count per
# institution are the confounders: T32-holding departments are, a priori,
# the larger and more research-intensive departments.
inst_base_funding = defaultdict(float)
inst_base_r01pi = defaultdict(set)
for _, row in dedup.iterrows():
    fy = row["Fiscal Year"]
    if pd.isna(fy):
        continue
    fy = int(fy)
    if 2010 <= fy <= 2020:
        inst_base_funding[row["Organization Name"]] += row["cost"]
        if row["Activity"] == "R01":
            inst_base_r01pi[row["Organization Name"]].add(row["Contact PI Person ID"])

reg_rows = []
for o in all_orgs:
    reg_rows.append({
        "org": o,
        "r01pi_late": len(inst_r01_pis_late.get(o, set())),
        "t32": 1 if o in inst_t32_years else 0,
        "log_base_funding": math.log1p(inst_base_funding.get(o, 0.0)),
        "base_r01pi": len(inst_base_r01pi.get(o, set())),
    })
reg_df = pd.DataFrame(reg_rows)
t32_model = {}
try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    # Unadjusted negative binomial
    m0 = smf.glm("r01pi_late ~ t32", data=reg_df,
                 family=sm.families.NegativeBinomial()).fit()
    # Adjusted for baseline funding and baseline R01-PI count
    m1 = smf.glm("r01pi_late ~ t32 + log_base_funding + base_r01pi", data=reg_df,
                 family=sm.families.NegativeBinomial()).fit()
    irr0 = math.exp(m0.params["t32"])
    irr1 = math.exp(m1.params["t32"])
    ci0 = m0.conf_int().loc["t32"].apply(math.exp)
    ci1 = m1.conf_int().loc["t32"].apply(math.exp)
    # Full adjusted-model coefficient table (all terms).
    ci1_all = m1.conf_int()
    term_label = {"Intercept": "Intercept",
                  "t32": "T32 grant held (FY2010-2020)",
                  "log_base_funding": "Log baseline FY2010-2020 NIH funding",
                  "base_r01pi": "Baseline R01-PI count (FY2010-2020)"}
    adj_terms = []
    for term in m1.params.index:
        adj_terms.append({
            "term": term_label.get(term, term),
            "coef": float(m1.params[term]),
            "exp_coef": float(math.exp(m1.params[term])),
            "ci_lo": float(math.exp(ci1_all.loc[term, 0])),
            "ci_hi": float(math.exp(ci1_all.loc[term, 1])),
            "p": float(m1.pvalues[term]),
        })
    # Poisson comparator -> Pearson dispersion statistic (justifies the
    # negative-binomial family over Poisson).
    mp = smf.glm("r01pi_late ~ t32 + log_base_funding + base_r01pi", data=reg_df,
                 family=sm.families.Poisson()).fit()
    pois_disp = float(mp.pearson_chi2 / mp.df_resid)
    t32_model = {
        "unadjusted_irr": irr0, "unadjusted_p": float(m0.pvalues["t32"]),
        "unadjusted_irr_ci": [float(ci0[0]), float(ci0[1])],
        "adjusted_irr": irr1, "adjusted_p": float(m1.pvalues["t32"]),
        "adjusted_irr_ci": [float(ci1[0]), float(ci1[1])],
        "adjusted_covariates": "log baseline FY2010-2020 funding + baseline R01-PI count",
        "adjusted_full_terms": adj_terms,
        "n_institutions": int(len(reg_df)),
        "adjusted_loglik": float(m1.llf),
        "adjusted_aic": float(m1.aic),
        "adjusted_deviance": float(m1.deviance),
        "adjusted_df_resid": int(m1.df_resid),
        "poisson_pearson_dispersion": pois_disp,
    }
    print(f"    Negative-binomial model, outcome = distinct R01 PIs FY2021-2025:")
    print(f"      Unadjusted T32 incidence rate ratio (IRR) = {irr0:.2f} "
          f"(p = {m0.pvalues['t32']:.4g})")
    print(f"      Adjusted T32 IRR = {irr1:.2f} "
          f"(95% CI {ci1[0]:.2f}-{ci1[1]:.2f}, p = {m1.pvalues['t32']:.4g}) "
          f"[adjusted for log baseline funding + baseline R01-PI count]")
except Exception as e:  # pragma: no cover
    print(f"    [statsmodels unavailable: {e}]")

# ----------------------------------------------------------------------
# 7. TIER A  --  INSTITUTIONAL SUSTAINABILITY / SURVIVAL
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("7. TIER A: INSTITUTIONAL FUNDING SUSTAINABILITY")
print("=" * 72)

inst_years = defaultdict(set)
for _, row in dedup.iterrows():
    org = row["Organization Name"]
    fy = row["Fiscal Year"]
    if pd.notna(org) and pd.notna(fy):
        inst_years[org].add(int(fy))

inst_first = {o: min(ys) for o, ys in inst_years.items()}

# Year-over-year retention: of departments funded in year Y, the fraction
# also funded in Y+1 (FY2015-2025).
yoy_retention = {}
for y in range(ANALYSIS_START, ANALYSIS_END):
    fy_set = {o for o, ys in inst_years.items() if y in ys}
    nxt = {o for o, ys in inst_years.items() if (y + 1) in ys}
    yoy_retention[y] = 100.0 * len(fy_set & nxt) / len(fy_set)
print(f"  Year-over-year department retention FY2015-FY2025: "
      f"mean {np.mean(list(yoy_retention.values())):.1f}% "
      f"(range {min(yoy_retention.values()):.1f}-{max(yoy_retention.values()):.1f}%)")

# How much of the FY2015->FY2024 count growth is genuinely new departments
# vs. re-entry of previously funded departments?
ever_before_2015 = {o for o, y in inst_first.items() if y < 2015}
print(f"  Institutions ever funded before FY2015: n={len(ever_before_2015)}")
print(f"  Of 87 institutions, first-ever funded year distribution: "
      f"<=2010 n={sum(1 for y in inst_first.values() if y<=2010)}, "
      f"2011-2014 n={sum(1 for y in inst_first.values() if 2011<=y<=2014)}, "
      f"2015-2021 n={sum(1 for y in inst_first.values() if 2015<=y<=2021)}, "
      f">=2022 n={sum(1 for y in inst_first.values() if y>=2022)}")

# Entrant cohorts: institutions whose first-EVER funded year is in 2015-2021,
# i.e. genuinely new entrants observable for >=4 years.
entrants = {o: y for o, y in inst_first.items() if 2015 <= y <= 2021}
print(f"  New entrant institutions (first funded FY2015-2021): n={len(entrants)}")

# Persistence: fraction of subsequent years funded; funded in FY2025
persist_fracs = []
funded_2025 = 0
for o, y0 in entrants.items():
    later = [y for y in range(y0 + 1, 2026)]
    if not later:
        continue
    frac = len([y for y in later if y in inst_years[o]]) / len(later)
    persist_fracs.append(frac)
    if 2025 in inst_years[o]:
        funded_2025 += 1
print(f"  Mean post-entry years funded: {100*np.mean(persist_fracs):.1f}%")
print(f"  Entrants still funded in FY2025: {funded_2025}/{len(entrants)} "
      f"({100*funded_2025/len(entrants):.1f}%)")

# Established institutions (first funded <= 2010) still funded FY2025
established = {o: y for o, y in inst_first.items() if y <= 2010}
est_2025 = sum(1 for o in established if 2025 in inst_years[o])
print(f"  Established institutions (first funded <=FY2010): n={len(established)}; "
      f"funded in FY2025: {est_2025} ({100*est_2025/len(established):.1f}%)")

# FY2024 -> FY2025 contraction: who dropped out?
funded_2024 = {o for o, ys in inst_years.items() if 2024 in ys}
funded_2025_set = {o for o, ys in inst_years.items() if 2025 in ys}
dropped = funded_2024 - funded_2025_set
added = funded_2025_set - funded_2024
print(f"\n  FY2024 funded departments: {len(funded_2024)}")
print(f"  FY2025 funded departments: {len(funded_2025_set)}")
print(f"  Dropped out FY2024->FY2025: n={len(dropped)}; "
      f"newly funded: n={len(added)}")
drop_tenure = [len([y for y in inst_years[o] if y <= 2024]) for o in dropped]
drop_funding24 = []
for o in dropped:
    drop_funding24.append(
        dedup[(dedup["Organization Name"] == o) &
              (dedup["Fiscal Year"] == 2024)]["cost"].sum())
single = 0
dropped_detail = []
if drop_tenure:
    print(f"  Dropped departments: median funded-years through FY2024 = "
          f"{np.median(drop_tenure):.0f}; median FY2024 funding = "
          f"${np.median(drop_funding24):,.0f}")
    single = sum(1 for o in dropped
                 if len([y for y in inst_years[o] if y <= 2024]) == 1)
    print(f"  Of dropped departments, {single}/{len(dropped)} had been funded "
          f"in only a single prior year")
    # Anonymised per-department detail (departments labelled A-E by descending
    # FY2024 funding) for the supplementary table.
    detail = sorted(zip(dropped,
                        [len([y for y in inst_years[o] if y <= 2024]) for o in dropped],
                        drop_funding24),
                    key=lambda r: -r[2])
    for i, (_, tenure, f24) in enumerate(detail):
        dropped_detail.append({"label": f"Department {chr(65+i)}",
                               "funded_years_through_fy2024": int(tenure),
                               "fy2024_funding": float(f24),
                               "fy2024_grants": int(((dedup['Organization Name'] ==
                                   detail[i][0]) & (dedup['Fiscal Year'] == 2024)).sum())})
    print("  Dropped-department detail (anonymised):")
    for d in dropped_detail:
        print(f"    {d['label']}: {d['funded_years_through_fy2024']} funded years, "
              f"FY2024 ${d['fy2024_funding']:,.0f} ({d['fy2024_grants']} grant-year[s])")

# ----------------------------------------------------------------------
# 8. TIER A  --  FORMAL TREND & CHANGEPOINT TESTS
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("8. TIER A: TREND & CHANGEPOINT STATISTICS (real USD)")
print("=" * 72)

years = list(range(ANALYSIS_START, ANALYSIS_END + 1))
real_series = np.array([fy_real[y] for y in years])
log_real = np.log(real_series)


def mann_kendall(x):
    """Two-sided Mann-Kendall trend test. Returns (S, Z, p, tau)."""
    n = len(x)
    s = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            s += np.sign(x[j] - x[i])
    var_s = n * (n - 1) * (2 * n + 5) / 18.0
    if s > 0:
        z = (s - 1) / math.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / math.sqrt(var_s)
    else:
        z = 0.0
    from math import erfc
    p = erfc(abs(z) / math.sqrt(2))
    tau = s / (0.5 * n * (n - 1))
    return s, z, p, tau


mk = mann_kendall(real_series)
print(f"  Mann-Kendall (FY2015-FY2025): S={mk[0]}, Z={mk[1]:.2f}, "
      f"p={mk[2]:.4g}, tau={mk[3]:.3f}")
mk23 = mann_kendall(np.array([fy_real[y] for y in range(2015, 2024)]))
print(f"  Mann-Kendall (FY2015-FY2023): S={mk23[0]}, Z={mk23[1]:.2f}, "
      f"p={mk23[2]:.4g}, tau={mk23[3]:.3f}")

try:
    from scipy import stats
    # log-linear regression FY2015-2025
    sl, ic, r, pp, se = stats.linregress(years, log_real)
    print(f"  Log-linear OLS FY2015-FY2025: growth={100*(math.exp(sl)-1):.2f}%/yr, "
          f"R^2={r**2:.3f}, p={pp:.4g}")
    # FY2015-2023 trend, projection for FY2024-2025
    yr23 = list(range(2015, 2024))
    lr23 = np.log([fy_real[y] for y in yr23])
    sl2, ic2, r2, p2, se2 = stats.linregress(yr23, lr23)
    resid = lr23 - (ic2 + sl2 * np.array(yr23))
    s_err = np.sqrt(np.sum(resid ** 2) / (len(yr23) - 2))
    xbar = np.mean(yr23)
    sxx = np.sum((np.array(yr23) - xbar) ** 2)
    for y in (2024, 2025):
        pred = ic2 + sl2 * y
        pi = stats.t.ppf(0.975, len(yr23) - 2) * s_err * math.sqrt(
            1 + 1 / len(yr23) + (y - xbar) ** 2 / sxx)
        obs = log_real[years.index(y)]
        lo, hi = math.exp(pred - pi), math.exp(pred + pi)
        print(f"  FY{y}: observed ${math.exp(obs)/1e6:.1f}M vs FY2015-2023 trend "
              f"projection ${math.exp(pred)/1e6:.1f}M "
              f"(95% PI ${lo/1e6:.1f}-{hi/1e6:.1f}M) "
              f"-> {'WITHIN' if lo<=math.exp(obs)<=hi else 'OUTSIDE'} interval")
except Exception as e:  # pragma: no cover
    print(f"  [scipy unavailable: {e}]")

# Funded-department count trend
dept_counts = np.array([concentration[y]["n_departments"] for y in years])
mk_dept = mann_kendall(dept_counts)
print(f"  Mann-Kendall funded-department count: Z={mk_dept[1]:.2f}, "
      f"p={mk_dept[2]:.4g}, tau={mk_dept[3]:.3f}")

# ----------------------------------------------------------------------
# SAVE
# ----------------------------------------------------------------------
output = {
    "meta": {
        "source_file": os.path.basename(REPORTER_FILE),
        "reporter_accessed": "2026-03-08",
        "analysis_run": "2026-10-12",
        "cpi_base_year": CPI_BASE_YEAR,
        "raw_rows": int(n_raw),
        "deduplicated_records": int(n_dedup),
        "duplicate_rows_removed": int(dup_app),
        "unique_projects": int(dedup["Project Number"].nunique()),
        "unique_institutions": int(dedup["Organization Name"].nunique()),
        "unique_pis": int(dedup["Contact PI Person ID"].nunique()),
        "total_funding_raw": float(raw_total),
        "total_funding_dedup": float(dedup_total),
    },
    "fy_nominal": {str(k): v for k, v in sorted(fy_nominal.items())},
    "fy_real": {str(k): v for k, v in sorted(fy_real.items())},
    "fy_records": {str(k): v for k, v in sorted(fy_records.items())},
    "fy_projects": {str(k): v for k, v in sorted(fy_projects.items())},
    "cagr": {"nominal": nom_cagr, "real": real_cagr},
    "dedup_sensitivity": {
        "raw_nominal": {str(y): v["raw_nominal"] for y, v in dedup_sens.items()},
        "raw_real": {str(y): v["raw_real"] for y, v in dedup_sens.items()},
        "raw_real_cagr": raw_real_cagr, "raw_nominal_cagr": raw_nom_cagr,
        "raw_peak_year": int(raw_peak_year),
        "raw_peak_to_2025_pct": float(100 * (raw_fy_real[2025] / raw_fy_real[2023] - 1)),
        "dup_appids": int(n_dup_appids),
        "dup_exact_cost": int(exact_dups),
        "dup_costdiff": int(costdiff_dups),
        "dup_type_distribution": {str(k): int(v) for k, v in dup_type_dist.items()},
    },
    "multi_pi_pct": float(multi_pi_pct),
    "budget_share": {str(k): v for k, v in budget_share.items()},
    "nih_appropriation_m": NIH_APPROPRIATION_M,
    "nih_rpg_funding_usd": NIH_RPG_FUNDING_USD,
    "nih_rpg_real": {str(k): v for k, v in nih_rpg_real.items()},
    "nih_rpg_award_counts": NIH_RPG_AWARD_COUNTS,
    "nih_appropriation_real": {str(k): v for k, v in nih_appropriation_real.items()},
    "nih_benchmark": nih_benchmark,
    "concentration": {str(k): v for k, v in concentration.items()},
    "fy2025_lorenz_curve": {
        "department_share": [float(v) for v in fy2025_lorenz_departments],
        "funding_share": [float(v) for v in fy2025_lorenz_funding],
    },
    "top_institutions_cumulative_real": {
        org: float(val) for org, val in inst_cum.head(15).items()},
    "mechanism_mix": mech_mix,
    "r01_r37_share": float(100 * r01_r37 / win_total_real),
    "subspecialty_real_funding": {k: float(v) for k, v in
                                  sorted(sub_real.items(), key=lambda x: -x[1])},
    "subspecialty_counts": {k: int(v) for k, v in sub_count.items()},
    "subspecialty_fy_real": {s: {str(y): v for y, v in d.items()}
                             for s, d in sub_fy_real.items()},
    "subspecialty_real_funding_B": {k: float(v) for k, v in
                                    sorted(sub_real_B.items(), key=lambda x: -x[1])},
    "classifier_robustness": {
        "rank_spearman_A_vs_B": float(rho_rank) if rho_rank is not None else None,
        "mean_jaccard_A_vs_B": mean_jaccard,
        "exact_setmatch_pct_A_vs_B": float(exact_pct),
        "unclassified_pct_A": float(unclassified_A),
        "rank_order_A": order_A,
        "rank_order_B": order_B,
    },
    "pi_pipeline": {
        "cohort_n": n_cohort, "transitioned_n": n_trans,
        "transition_rate": float(100 * n_trans / n_cohort),
        "median_years_to_R": float(np.median(times)) if times else None,
        "iqr_years_to_R": [float(np.percentile(times, 25)),
                           float(np.percentile(times, 75))] if times else None,
        "wide_cohort_n": len(wide_cohort), "wide_transitioned_n": len(wide_trans),
        "wide_transition_rate": float(100 * len(wide_trans) / len(wide_cohort)),
        "cutoff_sensitivity": {str(k): v for k, v in pipeline_sensitivity.items()},
        "cutoff_rationale": ("FY2019 provides at least six fiscal years of "
                             "potential follow-up through FY2025, longer than "
                             "the observed median five-year transition time"),
        "mentored_K_cohort": len(k_cohort), "mentored_K_transitioned": len(k_trans),
        "mentored_K_rate": float(100 * len(k_trans) / max(1, len(k_cohort))),
        "fellowship_cohort": len(f_cohort), "fellowship_transitioned": len(f_trans),
        "fellowship_rate": float(100 * len(f_trans) / max(1, len(f_cohort))),
        "t32_institution_test": {
            "n_t32_inst": len(t32_orgs), "n_no_t32_inst": len(no_t32_orgs),
            "median_r01pi_t32": float(np.median(r01_t32)),
            "median_r01pi_no_t32": float(np.median(r01_no)),
            "mannwhitney_p": float(p_u) if p_u is not None else None,
            "spearman_rho": float(rho) if rho is not None else None,
            "spearman_p": float(p_rho) if p_rho is not None else None,
        },
        "t32_adjusted_model": t32_model,
        "pi_institution_mobility": pi_mobility,
        "top20_pi_type_audit": top_pi_audit,
        "top20_pi_type_counts": dict(top_pi_category_counts),
    },
    "sustainability": {
        "yoy_retention": {str(k): v for k, v in yoy_retention.items()},
        "yoy_retention_mean": float(np.mean(list(yoy_retention.values()))),
        "first_funded_le2010": sum(1 for y in inst_first.values() if y <= 2010),
        "first_funded_2011_2014": sum(1 for y in inst_first.values() if 2011 <= y <= 2014),
        "first_funded_2015_2021": sum(1 for y in inst_first.values() if 2015 <= y <= 2021),
        "first_funded_ge2022": sum(1 for y in inst_first.values() if y >= 2022),
        "n_entrants_2015_2021": len(entrants),
        "mean_postentry_years_funded_pct": float(100 * np.mean(persist_fracs)),
        "entrants_funded_2025": funded_2025,
        "entrants_funded_2025_pct": float(100 * funded_2025 / len(entrants)),
        "n_established": len(established),
        "established_funded_2025_pct": float(100 * est_2025 / len(established)),
        "funded_2024": len(funded_2024), "funded_2025": len(funded_2025_set),
        "dropped_2024_2025": len(dropped), "added_2024_2025": len(added),
        "dropped_median_tenure": float(np.median(drop_tenure)) if drop_tenure else None,
        "dropped_median_fy2024_funding": float(np.median(drop_funding24)) if drop_tenure else None,
        "dropped_single_year": int(single) if drop_tenure else None,
        "dropped_detail": dropped_detail,
    },
    "trend_tests": {
        "mann_kendall_2015_2025": {"S": mk[0], "Z": mk[1], "p": mk[2], "tau": mk[3]},
        "mann_kendall_2015_2023": {"S": mk23[0], "Z": mk23[1], "p": mk23[2], "tau": mk23[3]},
        "mann_kendall_dept_count": {"Z": mk_dept[1], "p": mk_dept[2], "tau": mk_dept[3]},
    },
}
with open(os.path.join(OUT_DIR, "analysis_output.json"), "w") as f:
    json.dump(output, f, indent=2, default=str)

# also write the deduplicated cohort for the reproducibility repo
dedup_export = dedup.drop(columns=["subspecialties", "cost"]).copy()
dedup_export.to_csv(os.path.join(OUT_DIR, "deduplicated_cohort.csv"), index=False)

print("\n" + "=" * 72)
print(f"SAVED: {os.path.join(OUT_DIR, 'analysis_output.json')}")
print(f"SAVED: {os.path.join(OUT_DIR, 'deduplicated_cohort.csv')}")
print("=" * 72)
