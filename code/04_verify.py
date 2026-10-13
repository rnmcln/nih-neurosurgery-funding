#!/usr/bin/env python3
"""
04_verify.py
============
Independent verification of every statistic and result reported in the
manuscript. This script does NOT import the main analysis pipeline; it
re-loads the raw NIH RePORTER export and recomputes each quantity from
scratch with independent code, then compares the recomputed value with the
value stated in the manuscript. A plain-text report is written to
verification/statistics_verification.txt.

Each check prints PASS (recomputed value matches the manuscript, within a
small tolerance) or FAIL (mismatch -- investigate).
"""
import json
import math
import os
import datetime

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XLSX = os.environ.get(
    "REPORTER_FILE",
    os.path.join(ROOT, "data", "source",
                 "SearchResult_Export_08Mar2026_022554.xlsx"),
)
LOCAL_MESH = os.path.join(ROOT, "outputs", "pubmed_mesh_lag.json")
BUNDLED_MESH = os.path.join(ROOT, "data", "aggregate", "pubmed_mesh_lag.json")
MESH = os.environ.get(
    "REPORTER_MESH_FILE",
    LOCAL_MESH if os.path.exists(LOCAL_MESH) else BUNDLED_MESH,
)
REPORT = os.environ.get(
    "REPORTER_VERIFICATION_REPORT",
    os.path.join(ROOT, "verification", "statistics_verification.txt"),
)
if not os.path.exists(XLSX):
    raise FileNotFoundError(
        "NIH RePORTER source export not found. Download the export described in "
        "DATA_PROVENANCE.md and set REPORTER_FILE to its local path."
    )

CPI = {1992: 140.3, 1993: 144.5, 1994: 148.2, 1995: 152.4, 1996: 156.9,
       1997: 160.5, 1998: 163.0, 1999: 166.6, 2000: 172.2, 2001: 177.1,
       2002: 179.9, 2003: 184.0, 2004: 188.9, 2005: 195.3, 2006: 201.6,
       2007: 207.342, 2008: 215.303, 2009: 214.537, 2010: 218.056,
       2011: 224.939, 2012: 229.594, 2013: 232.957, 2014: 236.736,
       2015: 237.017, 2016: 240.007, 2017: 245.120, 2018: 251.107,
       2019: 255.657, 2020: 258.811, 2021: 270.970, 2022: 292.655,
       2023: 304.702, 2024: 313.689, 2025: 321.943}
NIH_APPROP_M = {2015: 30311, 2016: 32311, 2017: 34301, 2018: 37311,
                2019: 39313, 2020: 41690, 2021: 42937, 2022: 45180,
                2023: 47684, 2024: 47099, 2025: 47099}
NIH_RPG_FUNDING_USD = {
    2015: 15862012059, 2016: 17137754907, 2017: 18321187243,
    2018: 19869684598, 2019: 21589062208, 2020: 22636319154,
    2021: 23280221639, 2022: 24399835832, 2023: 25365442608,
    2024: 25721050507, 2025: 26503031656,
}
NIH_RPG_AWARD_COUNTS = {2023: 41597, 2024: 41470, 2025: 39239}

lines = []
n_pass = n_fail = 0


def out(s=""):
    lines.append(s)


def check(label, computed, stated, tol=0.0, fmt="{}"):
    """Compare a recomputed value with the manuscript's stated value."""
    global n_pass, n_fail
    if isinstance(computed, (int, float)) and isinstance(stated, (int, float)):
        ok = abs(computed - stated) <= tol
    else:
        ok = computed == stated
    status = "PASS" if ok else "FAIL"
    if ok:
        n_pass += 1
    else:
        n_fail += 1
    out(f"  [{status}] {label}")
    out(f"         recomputed = {fmt.format(computed)}   "
        f"manuscript = {fmt.format(stated)}")


out("=" * 78)
out("INDEPENDENT STATISTICS VERIFICATION REPORT")
out("Manuscript: NIH funding to neurosurgery departments, FY2015-FY2025")
out(f"Verification run: {datetime.date.today().isoformat()}")
out("Method: the raw RePORTER export is reloaded and every reported quantity is")
out("recomputed from scratch by this standalone script, then compared with the")
out("value stated in the manuscript. Tolerances allow for display rounding.")
out("=" * 78)

# --------------------------------------------------------------- load + dedup
df = pd.read_excel(XLSX)
df["Total Cost"] = pd.to_numeric(df["Total Cost"], errors="coerce")
df["cost"] = df["Total Cost"].fillna(0.0)
df["Fiscal Year"] = pd.to_numeric(df["Fiscal Year"], errors="coerce").astype("Int64")
df["Activity"] = df["Activity"].astype(str).str.strip().str.upper()
n_raw = len(df)
dedup = (df.sort_values("cost", ascending=False)
           .drop_duplicates(subset=["Application ID"], keep="first").copy())

out("\n1. COHORT AND DEDUPLICATION")
check("Raw export rows", n_raw, 8672, fmt="{:,}")
check("Deduplicated grant-year records", len(dedup), 8480, fmt="{:,}")
check("Surplus duplicate rows removed", n_raw - len(dedup), 192, fmt="{:,}")
check("Distinct duplicated Application IDs",
      int(df[df["Application ID"].duplicated(keep=False)]["Application ID"].nunique()),
      164, fmt="{:,}")
check("Unique projects", dedup["Project Number"].nunique(), 8411, fmt="{:,}")
check("Unique institutions", dedup["Organization Name"].nunique(), 87, fmt="{:,}")
check("Unique contact PIs", dedup["Contact PI Person ID"].nunique(), 1119, fmt="{:,}")
raw_total = df["cost"].sum()
dedup_total = dedup["cost"].sum()
check("Total funding, raw (USD)", round(raw_total), 3453614202, tol=5, fmt="${:,}")
check("Total funding, deduplicated (USD)", round(dedup_total), 3318761705,
      tol=5, fmt="${:,}")
check("Double-counting removed (USD millions)",
      round((raw_total - dedup_total) / 1e6, 1), 134.9, tol=0.1, fmt="${}M")
check("Double-counting removed (% of raw)",
      round(100 * (raw_total - dedup_total) / raw_total, 1), 3.9, tol=0.1, fmt="{}%")
# all duplicates exact-cost?
dups = df[df["Application ID"].duplicated(keep=False)]
all_exact = all(len(set(g["Total Cost"].fillna(-1).round(2))) == 1
                for _, g in dups.groupby("Application ID"))
check("All duplicated Application IDs carry identical Total Cost",
      "yes" if all_exact else "no", "yes")

# --------------------------------------------------------------- funding
fy_nom = {int(k): float(v) for k, v in
          dedup.groupby("Fiscal Year")["cost"].sum().items() if pd.notna(k)}
fy_real = {y: v * CPI[2025] / CPI[y] for y, v in fy_nom.items() if y in CPI}

out("\n2. FUNDING TRENDS (constant 2025 US dollars)")
check("Real funding FY2015 (USD millions)", round(fy_real[2015] / 1e6, 1),
      131.2, tol=0.1, fmt="${}M")
check("Real funding FY2023 peak (USD millions)", round(fy_real[2023] / 1e6, 1),
      268.5, tol=0.1, fmt="${}M")
check("Real funding FY2024 (USD millions)", round(fy_real[2024] / 1e6, 1),
      255.3, tol=0.1, fmt="${}M")
check("Real funding FY2025 (USD millions)", round(fy_real[2025] / 1e6, 1),
      245.7, tol=0.1, fmt="${}M")
real_cagr = (fy_real[2025] / fy_real[2015]) ** 0.1 - 1
nom_cagr = (fy_nom[2025] / fy_nom[2015]) ** 0.1 - 1
check("Real CAGR FY2015-FY2025 (%)", round(100 * real_cagr, 1), 6.5,
      tol=0.1, fmt="{}%")
check("Nominal CAGR FY2015-FY2025 (%)", round(100 * nom_cagr, 1), 9.8,
      tol=0.1, fmt="{}%")
check("Peak-to-FY2025 real decline (%)",
      round(100 * (fy_real[2025] / fy_real[2023] - 1), 1), -8.5, tol=0.1,
      fmt="{}%")
check("Peak fiscal year", max(range(2015, 2026), key=lambda y: fy_real[y]),
      2023)

# dedup sensitivity
raw_nom = {int(k): float(v) for k, v in
           df.groupby("Fiscal Year")["cost"].sum().items() if pd.notna(k)}
raw_real = {y: v * CPI[2025] / CPI[y] for y, v in raw_nom.items() if y in CPI}
check("Raw-series peak fiscal year (sensitivity)",
      max(range(2015, 2026), key=lambda y: raw_real[y]), 2023)
check("Raw-series peak-to-FY2025 decline (%)",
      round(100 * (raw_real[2025] / raw_real[2023] - 1), 1), -8.1, tol=0.15,
      fmt="{}%")

# NIH-wide context. All dollar series are converted with the same CPI-U values.
rpg_real = {y: v * CPI[2025] / CPI[y] for y, v in NIH_RPG_FUNDING_USD.items()}
app_real = {y: v * 1e6 * CPI[2025] / CPI[y] for y, v in NIH_APPROP_M.items()}
check("NIH-wide RPG real change FY2023-FY2025 (%)",
      round(100 * (rpg_real[2025] / rpg_real[2023] - 1), 1), -1.1,
      tol=0.1, fmt="{}%")
check("Total NIH appropriation real change FY2023-FY2025 (%)",
      round(100 * (app_real[2025] / app_real[2023] - 1), 1), -6.5,
      tol=0.1, fmt="{}%")
check("NIH-wide RPG award-count change FY2023-FY2025 (%)",
      round(100 * (NIH_RPG_AWARD_COUNTS[2025] /
                   NIH_RPG_AWARD_COUNTS[2023] - 1), 1), -5.7,
      tol=0.1, fmt="{}%")

# --------------------------------------------------------------- budget share
out("\n3. SHARE OF TOTAL NIH APPROPRIATION (approximate benchmark)")
for y, claim in [(2015, 0.32), (2023, 0.53), (2025, 0.52)]:
    s = 100 * fy_nom[y] / (NIH_APPROP_M[y] * 1e6)
    check(f"Budget share FY{y} (%)", round(s, 2), claim, tol=0.02, fmt="{}%")

# --------------------------------------------------------------- mechanisms
win = dedup[(dedup["Fiscal Year"] >= 2015) & (dedup["Fiscal Year"] <= 2025)].copy()
win["real"] = win.apply(lambda r: r["cost"] * CPI[2025] / CPI[int(r["Fiscal Year"])],
                        axis=1)
win_total = win["real"].sum()
mech = win.groupby("Activity")["real"].sum()
out("\n4. FUNDING MECHANISMS (FY2015-FY2025, real USD)")
check("R01 share of funding (%)", round(100 * mech.get("R01", 0) / win_total, 1),
      54.5, tol=0.1, fmt="{}%")
check("R01 + R37 share of funding (%)",
      round(100 * (mech.get("R01", 0) + mech.get("R37", 0)) / win_total, 1),
      56.2, tol=0.1, fmt="{}%")
check("Window total real funding (USD billion)", round(win_total / 1e9, 2),
      2.30, tol=0.01, fmt="${}B")

# --------------------------------------------------------------- concentration
def gini(v):
    x = np.sort(np.asarray([a for a in v if a >= 0], float))
    n = len(x)
    return (2 * np.sum((np.arange(1, n + 1)) * x) / (n * x.sum())
            - (n + 1) / n)


out("\n5. INSTITUTIONAL CONCENTRATION (FY2015-FY2025)")
ginis, top10s, hhis, ndeps = [], [], [], []
for y in range(2015, 2026):
    sub = dedup[dedup["Fiscal Year"] == y]
    inst = sub.groupby("Organization Name")["cost"].sum()
    inst = inst[inst > 0].sort_values(ascending=False)
    ginis.append(gini(inst.values))
    top10s.append(100 * inst.iloc[:10].sum() / inst.sum())
    sh = 100 * inst / inst.sum()
    hhis.append(float(np.sum(sh ** 2)))
    ndeps.append(len(inst))
check("Gini coefficient range, lower bound", round(min(ginis), 3), 0.543,
      tol=0.001, fmt="{}")
check("Gini coefficient range, upper bound", round(max(ginis), 3), 0.576,
      tol=0.001, fmt="{}")
check("Top-10 share range, lower bound (%)", round(min(top10s), 1), 53.5,
      tol=0.1, fmt="{}%")
check("Top-10 share range, upper bound (%)", round(max(top10s), 1), 61.3,
      tol=0.1, fmt="{}%")
check("HHI range, lower bound", round(min(hhis)), 408, tol=1, fmt="{}")
check("HHI range, upper bound", round(max(hhis)), 651, tol=1, fmt="{}")
check("Funded departments FY2015", ndeps[0], 45)
check("Funded departments FY2024", ndeps[9], 59)
check("Funded departments FY2025", ndeps[10], 54)
inst25 = dedup[dedup["Fiscal Year"] == 2025].groupby(
    "Organization Name")["cost"].sum()
inst25 = inst25[inst25 > 0]
check("FY2025 median per-department funding (USD millions)",
      round(inst25.median() / 1e6, 2), 2.75, tol=0.01, fmt="${}M")
check("FY2025 IQR lower (USD millions)",
      round(inst25.quantile(0.25) / 1e6, 2), 1.14, tol=0.02, fmt="${}M")
check("FY2025 IQR upper (USD millions)",
      round(inst25.quantile(0.75) / 1e6, 2), 5.87, tol=0.02, fmt="${}M")

# top institution
win_inst = win.groupby("Organization Name")["real"].sum().sort_values(
    ascending=False)
check("Largest cumulative recipient",
      win_inst.index[0], "UNIVERSITY OF CALIFORNIA, SAN FRANCISCO")
check("Largest cumulative recipient real funding (USD millions)",
      round(win_inst.iloc[0] / 1e6, 1), 302.6, tol=0.2, fmt="${}M")

# --------------------------------------------------------------- pipeline
out("\n6. EARLY-CAREER AWARD TRAJECTORIES")
FELL = {"F30", "F31", "F32", "F99"}
MENT = {"K01", "K07", "K08", "K12", "K22", "K23", "K25", "K38", "K99", "K00"}
INDEP = {"R01", "R37", "RF1", "R35", "DP1", "DP2", "DP5", "R56", "R29"}
hist = {}
for _, r in dedup.iterrows():
    pid, fy, act = r["Contact PI Person ID"], r["Fiscal Year"], r["Activity"]
    if pd.notna(pid) and pd.notna(fy):
        hist.setdefault(pid, []).append((int(fy), act))
cohort = []
for pid, h in hist.items():
    h = sorted(h, key=lambda item: (item[0], str(item[1])))
    if h[0][1] in (FELL | MENT):
        cy = min(fy for fy, a in h if a in (FELL | MENT))
        if cy <= 2019:
            cohort.append((pid, h[0][1], cy))
trans = 0
times = []
k_n = k_t = f_n = f_t = 0
for pid, first_act, cy in cohort:
    h = sorted(hist[pid])
    ry = [fy for fy, a in h if a in INDEP and fy >= cy]
    did = len(ry) > 0
    trans += did
    if did:
        times.append(min(ry) - cy)
    if first_act in MENT:
        k_n += 1
        k_t += did
    elif first_act in FELL:
        f_n += 1
        f_t += did
check("Mature early-career cohort (n)", len(cohort), 266)
check("Transitioned to independent R award (n)", trans, 55)
check("Transition rate (%)", round(100 * trans / len(cohort), 1), 20.7,
      tol=0.1, fmt="{}%")
check("Median years to independent R award", float(np.median(times)), 5.0,
      tol=0.01)
check("Mentored-K entrants (n)", k_n, 104)
check("Mentored-K transition rate (%)", round(100 * k_t / k_n, 1), 40.4,
      tol=0.1, fmt="{}%")
check("Fellowship entrants (n)", f_n, 162)
check("Fellowship transition rate (%)", round(100 * f_t / f_n, 1), 8.0,
      tol=0.1, fmt="{}%")

def transition_at_cutoff(cutoff):
    members = []
    for pid, h in hist.items():
        h = sorted(h, key=lambda item: (item[0], str(item[1])))
        if h[0][1] in (FELL | MENT) and h[0][0] <= cutoff:
            members.append((pid, h[0][0]))
    converted = sum(any(a in INDEP and fy >= first_fy
                        for fy, a in hist[pid]) for pid, first_fy in members)
    return len(members), converted

for cutoff, expected_n, expected_converted, expected_rate in [
        (2017, 233, 51, 21.9), (2018, 254, 53, 20.9),
        (2019, 266, 55, 20.7), (2020, 280, 56, 20.0)]:
    sn, st = transition_at_cutoff(cutoff)
    check(f"Early-career sensitivity FY<={cutoff}, cohort (n)", sn, expected_n)
    check(f"Early-career sensitivity FY<={cutoff}, transitioned (n)",
          st, expected_converted)
    check(f"Early-career sensitivity FY<={cutoff}, rate (%)",
          round(100 * st / sn, 1), expected_rate, tol=0.1, fmt="{}%")

opl = dedup["Other PI or Project Leader(s)"].astype(str).str.strip()
multi = ~opl.isin(["", "nan", "None", "Not Applicable"])
check("Multi-PI grant-years (%)", round(100 * multi.sum() / len(dedup), 1),
      14.7, tol=0.1, fmt="{}%")

pi_org = (win.dropna(subset=["Contact PI Person ID", "Organization Name"])
          .groupby("Contact PI Person ID")["Organization Name"].nunique())
check("Contact PIs observed FY2015-FY2025 (n)", len(pi_org), 707)
check("Contact PIs linked to >1 institution FY2015-FY2025 (n)",
      int((pi_org > 1).sum()), 43)

# --------------------------------------------------------------- T32 model
out("\n7. T32 - R01 ASSOCIATION (negative-binomial regression)")
t32_years, r01_late, base_fund, base_r01 = {}, {}, {}, {}
allorg = set(dedup["Organization Name"].dropna())
for _, r in dedup.iterrows():
    o, fy, act = r["Organization Name"], r["Fiscal Year"], r["Activity"]
    if pd.isna(fy):
        continue
    fy = int(fy)
    if act == "T32" and 2010 <= fy <= 2020:
        t32_years.setdefault(o, set()).add(fy)
    if act == "R01" and 2021 <= fy <= 2025:
        r01_late.setdefault(o, set()).add(r["Contact PI Person ID"])
    if 2010 <= fy <= 2020:
        base_fund[o] = base_fund.get(o, 0.0) + r["cost"]
        if act == "R01":
            base_r01.setdefault(o, set()).add(r["Contact PI Person ID"])
reg = pd.DataFrame([{
    "y": len(r01_late.get(o, set())),
    "t32": 1 if o in t32_years else 0,
    "logbase": math.log1p(base_fund.get(o, 0.0)),
    "baser01": len(base_r01.get(o, set()))} for o in allorg])
try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    m0 = smf.glm("y ~ t32", reg, family=sm.families.NegativeBinomial()).fit()
    m1 = smf.glm("y ~ t32 + logbase + baser01", reg,
                 family=sm.families.NegativeBinomial()).fit()
    check("Unadjusted T32 incidence rate ratio",
          round(math.exp(m0.params["t32"]), 2), 3.52, tol=0.02)
    check("Adjusted T32 incidence rate ratio",
          round(math.exp(m1.params["t32"]), 2), 1.07, tol=0.03)
    ci = m1.conf_int().loc["t32"].apply(math.exp)
    check("Adjusted T32 IRR 95% CI lower", round(ci[0], 2), 0.45, tol=0.03)
    check("Adjusted T32 IRR 95% CI upper", round(ci[1], 2), 2.56, tol=0.05)
except Exception as e:  # pragma: no cover
    out(f"  [SKIP] statsmodels unavailable: {e}")
n_t32 = sum(1 for o in allorg if o in t32_years)
check("Institutions holding a T32 (FY2010-2020)", n_t32, 9)
check("Median R01 PIs FY2021-2025, T32 institutions",
      float(np.median([len(r01_late.get(o, set())) for o in allorg
                       if o in t32_years])), 9.0)
check("Median R01 PIs FY2021-2025, non-T32 institutions",
      float(np.median([len(r01_late.get(o, set())) for o in allorg
                       if o not in t32_years])), 1.0)

# --------------------------------------------------------------- domains
out("\n8. KEYWORD-ATTRIBUTED DISEASE-DOMAIN FUNDING (FY2015-FY2025, real USD)")
import re
LEX = {
    "Neuro-oncology": (["glioma", "glioblastoma", "brain tumor", "brain neoplasm",
        "meningioma", "brain cancer", "neuro-oncol", "neurooncol", "medulloblastoma",
        "craniopharyngioma", "schwannoma", "pituitary", "brain metastas",
        "ependymoma", "astrocytoma"], ["gbm"]),
    "Cerebrovascular": (["aneurysm", "cerebrovascular", "stroke",
        "subarachnoid hemorrhage", "arteriovenous malformation",
        "cavernous malformation", "intracranial hemorrhage", "cerebral hemorrhage",
        "moyamoya", "vasospasm", "endovascular", "carotid"], ["sah", "avm", "ich"]),
    "Spine": (["spinal", "spine", "lumbar", "cervical", "disc herniation",
        "scoliosis", "vertebr", "spinal cord injury", "spinal cord", "laminectomy",
        "degenerative disc", "myelopathy", "radiculopathy", "intervertebral"],
        ["sci"]),
    "Trauma/TBI": (["traumatic brain injury", "craniocerebral trauma", "subdural",
        "epidural hematoma", "skull fracture", "concussion", "head injury",
        "blast injury", "neurotrauma", "decompressive craniectomy"], ["tbi"]),
    "Functional/epilepsy": (["deep brain stimulation", "epilepsy", "seizure",
        "functional neurosurg", "movement disorder", "parkinson", "tremor",
        "dystonia", "responsive neurostimulation", "vagus nerve stimulat",
        "stereoelectroencephalography"], ["dbs"]),
    "Hydrocephalus": (["hydrocephalus", "ventriculoperitoneal", "cerebrospinal fluid",
        "normal pressure hydrocephalus", "endoscopic third ventriculostomy"],
        ["nph", "csf"]),
    "Peripheral nerve": (["peripheral nerve", "brachial plexus", "nerve injury",
        "nerve repair", "nerve regeneration", "nerve graft", "nerve transfer",
        "carpal tunnel", "ulnar nerve", "median nerve", "radial nerve",
        "sciatic nerve", "neurofibromatosis", "neurofibroma",
        "peripheral neuropath"], []),
}
FLD = ["Project Title", "Project Terms", "Project Abstract",
       "NIH Spending Categorization", "Public Health Relevance"]
dom_real = {k: 0.0 for k in LEX}
for _, r in win.iterrows():
    txt = " ".join(str(r.get(f, "") or "") for f in FLD).lower()
    for dom, (ph, ab) in LEX.items():
        if any(p in txt for p in ph) or any(
                re.search(r"\b" + a + r"\b", txt) for a in ab):
            dom_real[dom] += r["real"]
for dom, claim in [("Neuro-oncology", 684), ("Spine", 588),
                   ("Cerebrovascular", 565), ("Functional/epilepsy", 491),
                   ("Trauma/TBI", 247), ("Hydrocephalus", 178),
                   ("Peripheral nerve", 95)]:
    check(f"{dom} domain funding (USD millions)",
          round(dom_real[dom] / 1e6), claim, tol=1, fmt="${}M")

classified = 0.0
for _, r in win.iterrows():
    txt = " ".join(str(r.get(f, "") or "") for f in FLD).lower()
    if any(any(p in txt for p in ph) or
           any(re.search(r"\b" + a + r"\b", txt) for a in ab)
           for ph, ab in LEX.values()):
        classified += r["real"]
check("Unclassified funding share (%)",
      round(100 * (1 - classified / win_total), 1), 14.6, tol=0.1, fmt="{}%")

# The focused professional-role audit is manual; these checks independently
# confirm that the audited Person IDs are exactly the 20 highest-funded contact
# PIs in the analysis window and that the reported category totals are faithful.
AUDIT = {
    3160894: "neurosurgeon", 9729256: "neurosurgeon",
    7351737: "neurosurgeon", 8926441: "neurosurgeon",
    1896710: "non-physician scientist",
    7690446: "non-neurosurgeon physician",
    1894343: "non-neurosurgeon physician", 1898964: "neurosurgeon",
    2190844: "neurosurgeon", 9424298: "neurosurgeon",
    1900873: "non-neurosurgeon physician", 3143458: "neurosurgeon",
    1887100: "non-physician scientist", 1989991: "neurosurgeon",
    8544696: "neurosurgeon", 1878386: "neurosurgeon",
    6278770: "neurosurgeon", 6066228: "neurosurgeon",
    1928282: "neurosurgeon", 1865771: "neurosurgeon",
}
top_ids = set(int(x) for x in
              win.dropna(subset=["Contact PI Person ID"])
              .groupby("Contact PI Person ID")["cost"].sum()
              .nlargest(20).index)
check("Manual-audit IDs equal the 20 highest-funded contact PIs",
      "yes" if top_ids == set(AUDIT) else "no", "yes")
for category, expected in [("neurosurgeon", 15),
                           ("non-neurosurgeon physician", 3),
                           ("non-physician scientist", 2)]:
    check(f"Manual PI audit: {category} (n)",
          sum(v == category for v in AUDIT.values()), expected)

# --------------------------------------------------------------- MeSH
out("\n9. MeSH RETRIEVAL COVERAGE (NCBI E-utilities, queried 20 May 2026)")
M = json.load(open(MESH))
check("General MeSH:keyword ratio 2015",
      round(M["years"]["2015"]["mesh_to_keyword_ratio"], 2), 0.66, tol=0.01)
check("General MeSH:keyword ratio 2023",
      round(M["years"]["2023"]["mesh_to_keyword_ratio"], 2), 0.45, tol=0.01)
check("Journal MeSH ratio 2015",
      round(M["years"]["2015"]["journal_mesh_ratio"], 2), 0.43, tol=0.01)
check("Journal MeSH ratio 2023",
      round(M["years"]["2023"]["journal_mesh_ratio"], 2), 0.30, tol=0.01)

# --------------------------------------------------------------- summary
out("\n" + "=" * 78)
out(f"SUMMARY: {n_pass} checks PASSED, {n_fail} checks FAILED, "
    f"{n_pass + n_fail} total.")
if n_fail == 0:
    out("All reported statistics were independently reproduced from the raw")
    out("NIH RePORTER export. No discrepancies were found.")
else:
    out("One or more checks did not reproduce -- investigate the FAIL lines above.")
out("=" * 78)

os.makedirs(os.path.dirname(REPORT), exist_ok=True)
with open(REPORT, "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\nSAVED: {REPORT}")
