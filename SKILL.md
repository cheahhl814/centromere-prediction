---
name: nf-centromere-agents
description: 'Run and interpret centromere annotation and analysis with the nf-centromere-agents pipeline (BettaMt-convention meta-skill): TideHunter monomer discovery from reads, TRF tandem-repeat density on the assembly, stringDecomposer monomer decomposition, per-scaffold summary. Use when the user wants to find/annotate candidate centromere regions on an ONT fish assembly, run the centromere Nextflow workflow, or interpret existing results/centromere/ outputs with healthy-range targets and red-flag diagnostics. Triggers: "centromere", "centromere analysis", "find centromeres", "centromere annotation", "TideHunter", "TRF density", "stringdecomposer", "centromere candidates", "betta centromere", "check centromere results".'
version: 1.0.0
updated: "2026-10-21"
triggers:
  - "centromere"
  - "centromere analysis"
  - "find centromeres"
  - "centromere annotation"
  - "TideHunter"
  - "TRF density"
  - "stringdecomposer"
  - "centromere candidates"
  - "betta centromere"
  - "check centromere results"
---

# nf-centromere-agents — Centromere Annotation & Analysis

> **v1.0.0.** A BettaMt-convention meta-skill for centromere annotation of ONT-based fish genome assemblies. Wraps the standalone Nextflow DSL2 pipeline `runners/nextflow-runner/` (5 processes: `tidehunter_cluster`, `trf_density`, `identify_regions`, `stringdecomposer`, `centromere_summary`) and teaches the agent how to run it and — more importantly — how to *interpret* the outputs with evidence-based verdicts.

**History:** the centromere module started as `modules/centromere.nf` inside Betta_v3, where it was an orphan (not wired into `main.nf`). It was separated into its own standalone pipeline `nf-centromere-agents` and **5 latent bugs were fixed** en route (see §4 below; the 5th — a Nextflow channel-semantics bug that silently limited `stringdecomposer` to a single execution — surfaced on the first real PR001 (*Parosphromenus rubrimontis*) run and was fixed in commit `143a5dc`). This skill is the BettaMt-format rebuild of that pipeline: the Nextflow code is kept verbatim, only the documentation layer was restructured.

## Audience

This skill serves two purposes:

- **AI Agents**: Triggered by the phrases above. The agent must follow the strict evidence chain (Preflight → Run → QC → Debug) and the Go/No-Go gates.
- **Human Users**: Provides a transparent, reproducible pipeline with interpretation guidance.

## When to Use This Skill

Use this skill when you need to:

- Find / annotate candidate centromere regions on a fish assembly (scaffolded or polished FASTA).
- Run the centromere workflow on trimmed nuclear reads + an assembly after a `Betta_v3`-style WGS run.
- Interpret existing `results/centromere/` outputs (monomers, density windows, decompositions, summary TSV).

**Do NOT use this skill** if:

- The organism is holocentric (e.g. *C. elegans*) or its centromere monomers are far outside 50–250 bp (e.g. *Drosophila* ~5–10 bp) — the hardcoded TRF filter is vertebrate-tuned (see §10 of `qc/results-qc`).
- You are annotating structural (gene) features, not centromeres — use a genome-annotation skill.
- You only want to *re-run* a pipeline whose `results/centromere/` is already complete and interpreted — just read the report.

## 0. Orchestrator — detect stage, route to the right sub-skill

This meta-skill is a **router**, not a doer. Its job is to ask: **"what stage is the user at, and which sub-skill should they invoke next?"**

### 0.1 Locate the run directory

The run directory is where `results/centromere/`, `params.json`, and `centromere_run.log` live — *not* the skill directory. Default: the current working directory. Override with `RUN_DIR` (analysis run dir; distinct from the skill's install path).

### 0.2 Detect the user's stage

```bash
# Stage detection ladder — first match wins (run in the analysis RUN_DIR)
test -f "$RUN_DIR/centromere_report.md"          && STAGE="qc-done"
test -f "$RUN_DIR/results/centromere/centromere_summary.tsv" && STAGE="qc"
test -f "$RUN_DIR/.nextflow" && test -f "$RUN_DIR/centromere_run.log" && STAGE="run-running"
test -f "$RUN_DIR/params.json" && test -f "$RUN_DIR/preflight_passed.md" && STAGE="run"
: "${STAGE:=preflight}"
```

If auto-detection is ambiguous, ask the user one short question (SP0).

### 0.5 Master ask-user stop point

**SP0 — Entry-stage ambiguity** (Evidence + Recommend + Options): if the stage ladder returns `preflight` but a `results/centromere/` directory already exists (stale results from a previous run), warn the user: `publishDir` uses `overwrite: false`, so stale files are silently kept. Recommend moving old results to trash before running. Confirm before proceeding.

### 0.3 Route to the right sub-skill

| Stage | Action |
| --- | --- |
| `preflight` | Invoke `preflight/run-preflight`. Validates reads/assembly/tool env, writes `params.json`. |
| `run` | Invoke `run/pipeline-run`. Runs the Nextflow workflow (or the thin bash recipe). |
| `qc` | Invoke `qc/results-qc`. Interprets outputs, writes `centromere_report.md`. |
| `run` failed | Invoke `debug/run-debug`. Matches the failure against the signature library. |
| `qc-done` | Surface `centromere_report.md` to the user. |

**Do not skip preflight.** Wrong inputs (mtDNA-polluted reads, missing assembly) waste ~30 min of a 24-cpu run before failing downstream.

## 0.4 Execution paths

| Path | When | Where |
|------|------|-------|
| **Nextflow runner** (default) | Every run — the pipeline *is* Nextflow; 5 published-output processes, `-resume`, trace/report/timeline | `runners/nextflow-runner/` |
| **Bash recipes** | Debugging a single stage manually (e.g. re-running TRF standalone) | Each sub-skill's `Procedure` section |

Both paths share one pixi manifest at the skill root (`pixi.toml`; run tools with `pixi run --manifest-path <skill-root>/pixi.toml <tool>`). Nextflow itself lives in the isolated `nextflow` env: `pixi run -e nextflow nextflow ...` (no-default-feature; only materialized on demand). The bash recipes remain the source of truth; the runner is a thin executor.

## A. Pipeline Architecture

The analysis is divided into **three sequential phases**: Preflight → Run → QC. The agent MUST complete each phase in order and pass the associated "Go/No-Go" gate.

### Phase 1: Preflight (The Audit)

**Goal**: Validate reads (nuclear/trimmed preferred — mtDNA pollutes the monomer set), validate the assembly, verify the tool env resolves, check disk, and write the run `params.json`.

- **Sub-skill**: `preflight/run-preflight`
- **Must-Verify**: inputs exist and are non-empty; pixi env solves; helper script present.
- **Exit Gate**: preflight report with verdict ≥ `GO-WITH-WARNINGS` + run `params.json`.

### Phase 2: Run (The Workflow)

**Goal**: Execute the 5-process pipeline. TideHunter (reads → tandem-repeat monomers, clustered via cd-hit-est at 80%) and TRF (assembly → 50–250 bp repeats → 100 kb density windows) run in parallel, then `identify_regions` picks the top-10 density windows per scaffold, `stringdecomposer` decomposes each candidate region against the monomers, and `centromere_summary` aggregates.

- **Sub-skill**: `run/pipeline-run`
- **Must-Verify**: preflight verdict ≥ `GO-WITH-WARNINGS`.
- **Exit Gate**: `results/centromere/centromere_summary.tsv` non-empty (or a documented no-signal verdict).

### Phase 3: QC (The Report)

**Goal**: Interpret every output with healthy-range targets and red flags, and write `centromere_report.md` with a verdict: `good` | `fragmented` | `monomer_filter_too_strict` | `no_centromeric_signal`.

- **Sub-skill**: `qc/results-qc`
- **Must-Verify**: run outputs present.
- **Exit Gate**: `centromere_report.md` exists with verdict.

## 4. The 5 bugs fixed in the extracted pipeline (read before trusting old results)

1. **Publish pattern typo** — was `pattern: 'consensus_clustered.fa'` for a file named `consensus_clustered.fasta` → nothing published.
2. **`candidate_regions.tsv` never published** — `identify_regions` had no `publishDir`.
3. **Bare `bedtools` in `centromere_summary`** — called outside the tool env; now wrapped in `pixi run`.
4. **Unused `val clust_id` input on `tidehunter_cluster`** — declared 2 inputs, workflow passed 1; would have crashed the first real run.
5. **`stringdecomposer` ran only once (queue-channel bug)** — `assembly` and `monomers` were single-value queue channels; Nextflow limits a multi-input process to the shortest channel. Fixed: `assembly` is a value channel (`Channel.value(file(...))`) and monomers are combined into the regions tuple (`.combine`). Verified on PR001: all 173 candidate regions run.

If the user still has old Betta_v3-era results, `debug/run-debug` covers the legacy quirks in its signature library.

## B. The Evidence Chain & Handoff Contract

| From → To | Artifact | Owner | Consumer |
| --- | --- | --- | --- |
| Upstream WGS run → Preflight | trimmed reads + polished assembly paths | user / upstream skill | `preflight/run-preflight` |
| Preflight → Run | `params.json` + `preflight_passed.md` | `preflight/run-preflight` | `run/pipeline-run` |
| Run → QC | `results/centromere/*` (6 artifact kinds) | `run/pipeline-run` | `qc/results-qc` |
| QC → User | `centromere_report.md` (verdict + warnings) | `qc/results-qc` | User |
| any failure → Debug | failing command's stderr | any phase | `debug/run-debug` |

The boundary between phases is the **filesystem**, not the agent's memory.

## C. Output contract — no artifact

This orchestrator produces no file of its own. All artifacts (`params.json`, `preflight_passed.md`, `centromere_report.md`) are produced by the sub-skills in the **analysis run dir**.

## D. Go/No-Go Gates

The agent must stop and warn the user if:

- **Preflight `NO-GO`**: inputs missing, tool env unsolvable, or < 5 GB disk.
- **Run silent no-signal**: `centromere_summary.tsv` empty is a *verdict* (`no_centromeric_signal` / `monomer_filter_too_strict`), not a crash — route to QC interpretation, never to "retry harder".
- **QC `NO-GO`**: `centromere_report.md` reports `NO-GO` (e.g. assembly fragmentation makes candidates uninterpretable).

## E. Common follow-ups

| User says | What to do |
| --- | --- |
| "Find centromeres for <species>" | Full chain: preflight → confirm inputs → run → QC → `centromere_report.md`. |
| "Check my centromere results" | Route to `qc/results-qc` with the existing `results/centromere/`. |
| "The run failed" | Route to `debug/run-debug` with the failing stderr / log tail. |
| "Widen the period range" | See `qc/results-qc` §Period range customization (post-process, don't fork the module). |
| "Is this skill up-to-date?" | `pixi run update-check` (see below). |

## Update Check

This skill ships with a self-update check that compares the deployed skill against the upstream `github.com/cheahhl814/nf-centromere-agents` repo via `git fetch` + SHA diff — no GitHub API call, no extra dependencies.

```bash
pixi run update-check
# Verdict legend (exit code in parens):
#   UP-TO-DATE       (0)  local HEAD matches origin/HEAD
#   LOCAL-AHEAD      (0)  unpushed local commits; no action needed
#   BEHIND-BY-N      (1)  upstream is N commits ahead → rsync from the source checkout
#   OFFLINE          (2)  git fetch failed (no network / no credentials); informational
#   NO-ORIGIN        (2)  no `origin` remote configured; informational
```

When `BEHIND-BY-N`, the script prints the canonical fix (rsync from the `@skills/<skill>/` source checkout to the deployed copy, then `diff -rq` to verify). When `OFFLINE`, it still prints `local_sha` + deployed `version:` so you can compare by hand. Full implementation: `bin/skill-update-check.py`.

## Glossary

- **Monomer**: consensus tandem repeat unit (50–250 bp target range) discovered from reads and clustered with cd-hit-est.
- **HOR**: Higher-Order Repeat — hierarchical monomer arrangement typical of vertebrate centromeres; visible as periodic dominant-monomer structure in the decomposition.
- **Density window**: 100 kb interval scored by fraction of bases covered by candidate repeats (column 7 of `centromere_density_100kb.bed`).
- **Candidate region**: top-10 density windows per scaffold merged into one span (`candidate_regions.tsv`) — a *candidate*, not a proven centromere.
- **`publishDir overwrite: false`**: stale outputs from earlier runs are silently kept; always move old `results/` aside before re-running.

## Verification

- [ ] Preflight performed; verdict ≥ `GO-WITH-WARNINGS`.
- [ ] Run completed; `results/centromere/centromere_summary.tsv` present.
- [ ] QC performed; `centromere_report.md` exists with verdict.

## Handoff pointers

After routing the user to a phase, state explicitly which sub-skill was invoked and what file it produces.

## Provenance

- Rebuilt with the [skill-forge](https://github.com/cheahhl814/skill-forge) meta-skill (`nf-centromere-agents.zip` as source; BettaMt pattern), pattern adopted from [BettaMt-agents](https://github.com/cheahhl814/BettaMt-agents), [Betta-WGS-agent](https://github.com/cheahhl814/Betta-WGS-agent), and [bacterial-genome-analysis](https://github.com/cheahhl814/bacterial-genome-analysis).
- Worked example: PR001 *Parosphromenus rubrimontis* — 500 Mb assembly, 173 contigs; TideHunter ~13 min + TRF ~26 min in parallel; stringdecomposer 173 tasks; full pipeline < 1 h on 24 cpus / 96 GB.