# Workstream 2 — RNAnneal Free Tier Launch

Owner(s): Lukas (product), Nadia (infrastructure)

Definition of Done (from master plan): soft-launch to 10–20 labs; stable operations with limits, landing page, and repeatable pipeline.

## Outcomes / Deliverables

- MVP scope defined (limits/presets/queue) + written data policy
- Landing page copy + onboarding path
- Soft launch to friendly labs; initial job intake and support loop
- Observability: logs, failures, runtime metrics; operational playbook
- Throughput + cost controls adequate for small-scale sustained usage
- Clear list of “paid tier” requirements informed by early users

## Timeline (from master plan, expanded)

### January (MVP definition)
- Define MVP scope: sequence limit(s), presets, queue behavior, retention, guardrails
- Draft data policy (storage, retention, privacy expectations)
- Create landing page copy (what it does, who it’s for, limitations, contact)

### February (Soft launch)
- Soft launch to 10–20 friendly labs
- Start collecting real jobs; establish a feedback + triage loop
- Implement observability (logs, failures, runtime); add runbooks for common failure modes

### March–May (Scale & polish)
- Improve throughput and cost controls (quotas, batching, queue discipline)
- Basic UX polish driven by failure cases + support burden
- Identify paid tier requirements and success criteria for conversion

## Dependencies

- Cloud stability and demo readiness (Workstream 1 synergy)
- Clear product boundaries so support load is bounded (limits + UX clarity)
- Validation stories (Workstream 3) improve adoption and outreach for Workstream 4

## Risks & Mitigations

- Demo/ops brittleness: prefer reliability over features; keep fallback walkthrough video for outreach and mid-term contexts.
- Support burden: enforce hard limits early; write an ops playbook; automate common checks.

## Near-term Actions (Next 6 Weeks)

- Weeks 3–4: make MVP build decisions (limits/presets/retention); plan observability baseline.
- February: soft launch to friendly labs; begin collecting jobs; iterate based on runtime/failure telemetry.

