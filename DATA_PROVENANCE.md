# Data provenance

## NIH RePORTER cohort

- Source: [NIH RePORTER](https://reporter.nih.gov/)
- Access date: 8 March 2026
- Export filename: `SearchResult_Export_08Mar2026_022554.xlsx`
- Query: exact Department field value `NEUROSURGERY`
- Fiscal-year restriction at retrieval: none
- Activity-code restriction at retrieval: none
- Funding-institute restriction at retrieval: none
- Export rows: 8,672
- SHA-256: `c6e63a2fc6a81b9cea7acfd79104f0ce3a7ed62758940eacb11464f83fbc0027`

The source export is not redistributed here because it contains row-level
identifiers and applicant-supplied text fields. Reproduction requires a fresh
download from NIH RePORTER or the researcher's lawfully retained export. Place
the file outside version control and set `REPORTER_FILE` to its absolute path.

The analysis retains, for each Application ID, the row with the largest
reported Total Cost. In this export, 164 Application IDs appeared more than
once; removing 192 surplus rows produced 8,480 grant-year records. The locally
generated deduplicated file is ignored by Git and must not be committed.

## PubMed retrieval counts

- Source: NCBI E-utilities `esearch` count endpoint
- Query date used in the study: 20 May 2026
- Publication years: 2015–2025
- Complete search strings: `code/02_pubmed_mesh_lag.py`

Counts for 2024–2025 are retained for completeness but were excluded from
trend inference because indexing of recent records is incomplete.

## Inflation and contextual benchmark series

- CPI-U: US Bureau of Labor Statistics series `CUUR0000SA0`; annual averages,
  accessed 20 May 2026. The 2025 average uses the 11 published monthly values.
- NIH Research Project Grant obligations and award counts: NIH Data Book,
  Research Project Grants, FY2015–FY2025; accessed 12 October 2026.
- Total NIH appropriations: NIH Office of Budget and Congressional Research
  Service report R43341. This series is used only as a contextual benchmark.

All precomputed files under `data/aggregate/` contain aggregate counts,
funding summaries, model estimates, or query totals. They contain no grant-
level narratives or row-level investigator records.
