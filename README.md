# NIH funding to neurosurgery departments, FY2015–FY2025

Reproducibility materials for an analysis of National Institutes of Health
(NIH) funding administratively assigned to neurosurgery departments, early-
career award trajectories, institutional training infrastructure, and MeSH
retrieval coverage.

This release contains analysis code, non-identifiable aggregate results,
figures, query documentation, and an independent verification report. It does
not redistribute the NIH RePORTER row-level export, applicant-supplied
narratives, investigator-level records, or a derived row-level cohort.

## Main findings

- Inflation-adjusted departmental funding rose from $131.2 million in FY2015 to
  $268.5 million in FY2023, then fell to $245.7 million in FY2025.
- The 8.5% departmental decline after FY2023 was steeper than the 1.1% real
  decline in NIH-wide Research Project Grant obligations over the same period.
- Institutional concentration was stable: the Gini coefficient ranged from
  0.543 to 0.576 and the ten largest recipients held 53.5%–61.3% of funding.
- Later independent R-type funding followed 40.4% of mentored-K entries and
  8.0% of fellowship entries. The unadjusted association between institutional
  T32 status and later R01 capacity disappeared after adjustment for baseline
  research intensity (adjusted incidence rate ratio 1.07, 95% CI 0.45–2.56).
- MeSH-restricted retrieval captured less than half as many neurosurgical
  records as a keyword-augmented search by 2023.

All funding results describe departments as administrative units. They do not
identify the training, degree, or clinical specialty of every principal
investigator.

## Repository contents

```text
code/
  01_analysis.py          Main analysis and aggregate-result generation
  02_pubmed_mesh_lag.py   PubMed/MeSH retrieval counts
  03_figures.py           Main and supplementary figures
  04_verify.py            Independent recomputation of reported statistics

data/aggregate/
  analysis_output.json    Aggregate results used by the figure script
  pubmed_mesh_lag.json    Saved PubMed counts from 20 May 2026

figures/                  Publication figures in PNG and TIFF
verification/             Independent report: 91 checks passed, 0 failed
DATA_PROVENANCE.md        Source, query, checksum, and retrieval details
DATA_LICENSE.md           Rights and redistribution statement
CITATION.cff              Software/repository citation metadata
```

## Reproduce the analysis

Python 3.10 or later is required.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Download the NIH RePORTER export described in
[DATA_PROVENANCE.md](DATA_PROVENANCE.md), then provide its local path. The file
is intentionally not included in this repository.

```bash
export REPORTER_FILE="/absolute/path/to/SearchResult_Export_08Mar2026_022554.xlsx"
python code/01_analysis.py
python code/04_verify.py
python code/03_figures.py
```

Generated files are written to the ignored `outputs/` directory. The analysis
script creates a local deduplicated cohort solely for reproducibility; it is not
tracked or distributed. To regenerate the live PubMed counts:

```bash
python code/02_pubmed_mesh_lag.py
```

An NCBI API key is optional and can be supplied as `NCBI_API_KEY`. PubMed
counts can drift as records and indexing are updated; the bundled aggregate
file records the counts used in the study.

To regenerate figures directly from the bundled aggregate results, without the
row-level source export:

```bash
REPORTER_OUTPUT_DIR=data/aggregate python code/03_figures.py
```

## Verification

`code/04_verify.py` is independent of the main analysis script. It reloads the
user-supplied export and recomputes cohort, funding, concentration,
early-career, regression, disease-domain, and retrieval results. The included
report records 91 passed checks and no failures.

## Data and rights

The code is licensed under the MIT License. Aggregate outputs and author-created
figures are released under CC BY 4.0. Third-party row-level source data are not
licensed by this repository and are not redistributed. See
[DATA_LICENSE.md](DATA_LICENSE.md) and
[DATA_PROVENANCE.md](DATA_PROVENANCE.md) for details.
