# Task: eSHAPE Benchmark (Workstream 1 — MII Grant)

Timeframe: February → May (data expected Mar 15–Apr 15; benchmark scaling April–May)

## Goal

Complete the 96‑RNA benchmark using experimental eSHAPE data (EclipseBio), producing analysis-quality results ready to support the final MII report.

## Inputs

- Final RNA target list + batching plan (locked)
- EclipseBio contract/quote and sequence submission pipeline
- Data ingestion + QC plan for eSHAPE outputs
- Scoring pipeline that can be re-run reproducibly on updated data drops

## Deliverables

- Integrated dataset (raw → processed) with provenance and QC summary
- Per-target and aggregate benchmark metrics (versioned, reproducible)
- M1/M2 analysis artifacts (plots/tables) suitable for report inclusion

## Milestones / Steps

- February: submit sequences; confirm timelines; prepare ingestion + QC scripts before data arrives
- March–April: ingest EclipseBio data on arrival; run initial benchmark pass; identify anomalies
- April–May: run benchmark at scale; iterate scoring model where justified (track changes + ablations)

## Risks & Mitigations

- Data arrival uncertainty: front-load tooling so analysis begins immediately upon receipt.
- Pipeline churn: version and freeze scoring definitions per report milestone; track deltas explicitly.

