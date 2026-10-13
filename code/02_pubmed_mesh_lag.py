#!/usr/bin/env python3
"""
02_pubmed_mesh_lag.py
=====================
Re-verifies the MeSH-indexing-lag analysis by querying NCBI E-utilities
(esearch, rettype=count) directly. Three retrieval strategies:

  Strategy 1  MeSH-restricted     : ("Neurosurgery"[MeSH] OR
                                     "Neurosurgical Procedures"[MeSH])
  Strategy 2  Keyword-augmented   : Strategy 1 OR title/abstract keywords
  Strategy 3  Journal subset      : nine core neurosurgery journals,
                                    with and without a MeSH restriction

The MeSH-to-keyword ratio (Strategy 1 / Strategy 2) and the journal-subset
MeSH-coverage ratio quantify how much of the neurosurgical literature a
MeSH-only search misses, and whether that gap is widening over time.

Output: outputs/pubmed_mesh_lag.json
"""

import json
import os
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.environ.get("REPORTER_OUTPUT_DIR", os.path.join(ROOT, "outputs"))
os.makedirs(OUT_DIR, exist_ok=True)
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

MESH = '("Neurosurgery"[MeSH] OR "Neurosurgical Procedures"[MeSH])'
KEYWORD_AUG = (
    '("Neurosurgery"[MeSH] OR "Neurosurgical Procedures"[MeSH] OR '
    'neurosurg*[tiab] OR craniotomy[tiab] OR craniectomy[tiab] OR '
    'laminectomy[tiab] OR discectomy[tiab] OR "spinal fusion"[tiab] OR '
    '"deep brain stimulation"[tiab] OR "ventriculoperitoneal shunt"[tiab])'
)
JOURNALS = (
    '("J Neurosurg"[Journal] OR "Neurosurgery"[Journal] OR '
    '"World Neurosurg"[Journal] OR "Neurosurg Focus"[Journal] OR '
    '"Acta Neurochir (Wien)"[Journal] OR "J Neurosurg Pediatr"[Journal] OR '
    '"J Neurosurg Spine"[Journal] OR "Stereotact Funct Neurosurg"[Journal] OR '
    '"Neurosurg Rev"[Journal])'
)


def count(term):
    """Return the PubMed result count for a query term."""
    request_params = {
        "db": "pubmed", "term": term, "retmax": 0, "rettype": "count"}
    if NCBI_API_KEY:
        request_params["api_key"] = NCBI_API_KEY
    params = urllib.parse.urlencode(request_params)
    url = f"{ESEARCH}?{params}"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                xml = resp.read().decode()
            c = xml.split("<Count>")[1].split("</Count>")[0]
            time.sleep(0.4)            # be polite to NCBI (<3 req/s)
            return int(c)
        except Exception as e:          # pragma: no cover
            time.sleep(2 * (attempt + 1))
            last = e
    raise RuntimeError(f"esearch failed for: {term[:60]} ({last})")


results = {"queried": time.strftime("%Y-%m-%d"), "years": {}}
print(f"{'Year':>6} {'MeSH':>8} {'Keyword':>9} {'ratio':>7} "
      f"{'JrnlAll':>8} {'JrnlMeSH':>9} {'Jratio':>7}")
for yr in range(2015, 2026):
    yf = f'"{yr}"[pdat]'
    mesh = count(f"{MESH} AND {yf}")
    kw = count(f"{KEYWORD_AUG} AND {yf}")
    jall = count(f"{JOURNALS} AND {yf}")
    jmesh = count(f"{JOURNALS} AND {MESH} AND {yf}")
    ratio = mesh / kw if kw else None
    jratio = jmesh / jall if jall else None
    results["years"][str(yr)] = {
        "mesh_restricted": mesh, "keyword_augmented": kw,
        "mesh_to_keyword_ratio": ratio,
        "journal_all": jall, "journal_mesh": jmesh,
        "journal_mesh_ratio": jratio,
    }
    print(f"{yr:>6} {mesh:>8} {kw:>9} {ratio:>7.3f} "
          f"{jall:>8} {jmesh:>9} {jratio:>7.3f}")

# Summary: ratio change across the fully-indexed window 2015-2023
y15 = results["years"]["2015"]
y23 = results["years"]["2023"]
results["summary"] = {
    "general_ratio_2015": y15["mesh_to_keyword_ratio"],
    "general_ratio_2023": y23["mesh_to_keyword_ratio"],
    "journal_ratio_2015": y15["journal_mesh_ratio"],
    "journal_ratio_2023": y23["journal_mesh_ratio"],
    "note": ("2024-2025 shown for completeness but excluded from trend claims "
             "because MeSH indexing of recent publications is incomplete."),
}
print("\nGeneral PubMed MeSH:keyword ratio  2015 -> 2023: "
      f"{y15['mesh_to_keyword_ratio']:.3f} -> {y23['mesh_to_keyword_ratio']:.3f}")
print("Journal-subset MeSH coverage ratio 2015 -> 2023: "
      f"{y15['journal_mesh_ratio']:.3f} -> {y23['journal_mesh_ratio']:.3f}")

with open(os.path.join(OUT_DIR, "pubmed_mesh_lag.json"), "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSAVED: {os.path.join(OUT_DIR, 'pubmed_mesh_lag.json')}")
