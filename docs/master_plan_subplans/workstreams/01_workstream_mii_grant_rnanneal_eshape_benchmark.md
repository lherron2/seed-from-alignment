# Workstream 1 — MII Grant (RNAnneal + eSHAPE Benchmark)

Owner(s): Lukas (primary), Nadia (cloud infrastructure support)

Definition of Done (from master plan): mid-term demo in February; full 96‑RNA benchmark with eSHAPE data + final report by June 30.

## Outcomes / Deliverables

- Live cloud demo suitable for MII mid-term presentation (Feb 16–20)
- Benchmark plan: RNA target list, batching plan, scoring methodology, and risk register
- EclipseBio execution: quote/contract, sequence submission, data ingestion plan
- Benchmark execution: 96 RNAs scored; M1/M2 analysis complete
- Final MII report + deliverables submitted by June 30

## Timeline (from master plan, expanded)

### January (Weeks 1–4)
- Lock RNA target list + batching plan (96 RNAs; include rationale + coverage strategy)
- Initiate EclipseBio contract/quote via AWS ($4,000) and confirm operational timeline
- Stabilize cloud deployment + a “boring reliable” demo URL; prepare fallback (video/screenshots)

### February (Weeks 5–8)
- Mid-term delivery: live demo + benchmark plan (Feb 16–20)
- Submit sequences to EclipseBio (as soon as contract path is cleared)
- Document risks explicitly (data slip, demo reliability, throughput/cost)

### March (Weeks 9–12)
- Integrate EclipseBio data when it arrives (expected Mar 15–Apr 15 window)
- Expand benchmark execution; formalize API endpoints + documentation needed for repeatability

### April–May (Weeks 13–21)
- Run benchmark at scale; iterate scoring model as needed
- Complete M1/M2 analysis; draft report sections continuously (avoid end-loaded writing)

### June (Weeks 22–24)
- Final report polish + submission by June 30
- Package artifacts (data, plots, evaluation tables, reproducibility notes)

## Key Dependencies

- Cloud demo stability (also supports Workstream 2/4 credibility)
- EclipseBio data arrival timing (do not block mid-term on experimental data)
- Clear definition of benchmark “success” and scoring approach (consistent across runs)

## Risks & Pre-decided Responses (from master plan)

- If EclipseBio slips relative to mid-term: treat mid-term as “experiment underway”; ensure demo + benchmark plan + preliminary scoring (public SHAPE data) are strong; present clear timeline.
- If cloud demo is brittle: ship a clean walkthrough video + screenshots; optimize for reliability over features.

## Near-term Actions (Next 6 Weeks)

- Weeks 1–2: finalize EclipseBio meeting/quote; stabilize cloud URL; finalize RNA list + batching plan; define mid-term live demo scope.
- Weeks 3–4: draft mid-term slides; run demo end-to-end incl failure modes; implement pitch deck feedback (cross-workstream).
- Week 5: mid-term rehearsal + backup plan; finalize preprint checklist (cross-workstream); ensure benchmark plan is presentation-ready.
- Week 6: deliver mid-term; convert feedback into a March–June checklist.

