---
name: centromere-prediction
description: 'Run and interpret centromere prediction for ONT-based fish genome assemblies (BettaMt-convention meta-skill wrapping a 5-process Nextflow DSL2 pipeline): TideHunter monomer discovery from reads, TRF tandem-repeat density on the assembly, stringDecomposer monomer decomposition, per-scaffold summary. Use when the user wants to find/annotate candidate centromere regions on a fish assembly, run the centromere Nextflow workflow, or interpret existing results/centromere/ outputs with healthy-range targets and red-flag diagnostics. Outputs are candidate regions, not proven centromeres. Triggers: "centromere", "centromere analysis", "centromere prediction", "find centromeres", "centromere annotation", "TideHunter", "TRF density", "stringdecomposer", "centromere candidates", "betta centromere", "check centromere results".'
version: 1.1.0
updated: "2026-11-03"
triggers:
  - "centromere"
  - "centromere analysis"
  - "centromere prediction"
  - "find centromeres"
  - "centromere annotation"
  - "TideHunter"
  - "TRF density"
  - "stringdecomposer"
  - "centromere candidates"
  - "betta centromere"
  - "check centromere results"
skills:
  - "read-qc-trimming"
  - "chartz"
---

# centromere-prediction — Centromere Prediction

> **v1.1.0.** A BettaMt-convention meta-skill for centromere annotation of ONT-based fish genome assemblies. Wraps the standalone Nextflow DSL2 pipeline `runners/nextflow-runner/` (5 processes: `tidehunter_cluster`, `trf_density`, `identify_regions`, `stringdecomposer`, `centromere_summary`) and teaches the agent how to run it and — more importantly — how to *interpret* the outputs with evidence-based verdicts.

**History:** the module began as an orphan `modules/centromere.nf` in Betta_v3, was split into the standalone `nf-centromere-agents` pipeline (5 latent bugs fixed, §4) and rebuilt here in the BettaMt format. **v1.1.0** (first real *Betta coccina* run) made the period range a parameter, bounded `stringdecomposer` input, fixed an unpublished decomposition table, and added satellite mapping, a telomere cross-check and evidence figures (§4c).

## Audience

**AI agents** follow the evidence chain (Preflight → Run → QC → Debug) and the Go/No-Go gates; **human users** get a reproducible pipeline with interpretation guidance.

## When to Use This Skill

Use it to:

- Find / annotate candidate centromere regions on a fish assembly (scaffolded or polished FASTA).
- Run the centromere workflow on trimmed nuclear reads + an assembly after a `Betta_v3`-style WGS run.
- Interpret existing `results/centromere/` outputs.

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

The **Nextflow runner** (`runners/nextflow-runner/`) is the default: 5 published-output processes, `-resume`, trace/report/timeline. **Bash recipes** (each sub-skill's `Procedure`) are for debugging one stage, and the QC helpers in `bin/` run outside Nextflow. Both share the pixi manifest at the skill root (`pixi run --manifest-path <skill-root>/pixi.toml <tool>`); Nextflow lives in the isolated `nextflow` env (`pixi run -e nextflow nextflow ...`).

## A. Pipeline Architecture

The analysis is divided into **three sequential phases**: Preflight → Run → QC. The agent MUST complete each phase in order and pass the associated "Go/No-Go" gate.

### Phase 1: Preflight (The Audit)

**Goal**: Validate reads (nuclear/trimmed preferred — mtDNA pollutes the monomer set), validate the assembly, verify the tool env resolves, check disk, and write the run `params.json`.

- **Sub-skill**: `preflight/run-preflight`
- **Must-Verify**: inputs exist and are non-empty; pixi env solves; helper script present.
- **Exit Gate**: preflight report with verdict ≥ `GO-WITH-WARNINGS` + run `params.json`.

### Phase 2: Run (The Workflow)

**Goal**: Execute the 5-process pipeline. TideHunter (reads → tandem-repeat monomers, clustered via cd-hit-est at 80%) and TRF (assembly → repeats of `--period_min`–`--period_max` bp, default 50–250 → 100 kb density windows) run in parallel, then `identify_regions` picks the densest `--region_window` (default 1 Mb) per scaffold, TideHunter's clusters are capped to the `--max_monomers` (default 100) most abundant, `stringdecomposer` decomposes each candidate region against the monomers, and `centromere_summary` aggregates.

- **Sub-skill**: `run/pipeline-run`
- **Must-Verify**: preflight verdict ≥ `GO-WITH-WARNINGS`.
- **Exit Gate**: `results/centromere/centromere_summary.tsv` non-empty (or a documented no-signal verdict).

### Phase 3: QC (The Report)

**Goal**: Interpret every output with healthy-range targets and red flags, and write `centromere_report.md` with a verdict: `good` | `partial` | `fragmented` | `monomer_filter_too_strict` | `no_centromeric_signal`. A negative verdict from the default 50–250 bp run is provisional: rerun wide (`--period_max 2000`) first. Then map the dominant satellite and cross-check telomeres (steps 6–8).

- **Sub-skill**: `qc/results-qc`
- **Must-Verify**: run outputs present.
- **Exit Gate**: `centromere_report.md` exists with verdict.

## 4. The 7 bugs fixed in the pipeline (read before trusting old results)

1. **Publish pattern typo** — was `pattern: 'consensus_clustered.fa'` for a file named `consensus_clustered.fasta` → nothing published.
2. **`candidate_regions.tsv` never published** — `identify_regions` had no `publishDir`.
3. **Bare `bedtools` in `centromere_summary`** — called outside the tool env; now wrapped in `pixi run`.
4. **Unused `val clust_id` input on `tidehunter_cluster`** — declared 2 inputs, workflow passed 1; would have crashed the first real run.
5. **`stringdecomposer` ran only once (queue-channel bug)** — `assembly` and `monomers` were single-value queue channels; Nextflow limits a multi-input process to the shortest channel. Fixed: `assembly` is a value channel (`Channel.value(file(...))`) and monomers are combined into the regions tuple (`.combine`). Verified on PR001: all 173 candidate regions run.
6. **`decomposition.tsv` was only the log** — stringDecomposer's stdout was redirected into it, and the real `decomposition_dir/final_decomposition.tsv` was never published. Now copied to `decomposition.tsv`.
7. **`stringdecomposer` OOM** — the merged top-10 span was often nearly a whole scaffold (2.6 Mb) and was aligned to every cd-hit cluster (3,017): SIGKILL, ~10 GB per 1 Mb × 100 monomers. Fixed in v1.1.0: densest-window regions + capped monomers (§4c). Also: BED start is 0-based but samtools is 1-based (off by one), and an empty monomer set segfaults the aligner (now skipped).

If the user still has old Betta_v3-era results, `debug/run-debug` covers the legacy quirks in its signature library.

### 4b. Adaptations made during the rebuild (4 — all verified by the `-profile test -stub-run` canary)

1. **Manifest paths**: `${baseDir}/envs/centromere/pixi.toml` → `${projectDir}/../../pixi.toml`; `bash ${baseDir}/bin/...` → `${projectDir}/../../bin/...`.
2. **Strict-DSL publishDir** (Nextflow 24+): the `stringdecomposer` bare-string `publishDir` failed with `No such property: scaffold for class: ProcessDslV1`; wrapped in a closure.
3. **Strict-DSL main.nf**: the top-level params `if` (hard error in Nextflow ≥26) moved inside `workflow`; `test` profile registered.
4. **Pixi-env untaint**: tasks launched via `pixi run -e nextflow` inherit `PIXI_*` vars, so the nested `pixi run --manifest-path` resolved nextflow's env (`Tool: command not found`). Every pixi call is prefixed with `env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE`.

### 4c. v1.1.0 changes (verified: `-profile test -stub-run` and the new real-run `-profile test_sat` on a satellite fixture)

- **Parameters** (`--period_min/--period_max/--min_copies/--region_window/--max_monomers`): TRF's max period is its own 7th argument, so widening needs a rerun, not a post-filter. TRF caps at 2000.
- **`bin/rank_monomers.py`** ranks cd-hit clusters by size (`monomer_ranking.tsv`, `monomers_top.fasta`); `bin/extract_candidate_region.sh` takes a window argument (0 = legacy merged span).
- **QC helpers** (steps 6–8): `bin/decomp_stats.py`, `satellite_map.py`, `telomere_crosscheck.py`, `evidence_tracks.py`. The env now includes `blast` and `tidk`.

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
| "Widen the period range" / negative narrow result | Rerun with `--period_max 2000` (TRF itself caps at `period_max`); see `qc/results-qc` §Period range. |
| "Which chromosomes carry the satellite?" / "centromere vs telomere" | `qc/results-qc` steps 6–8: `decomp_stats` → `satellite_map` → `telomere_crosscheck` → `evidence_tracks`. |
| "Is this skill up-to-date?" | `pixi run update-check` (see below). |

## Update Check

`pixi run update-check` compares the deployed skill with `github.com/cheahhl814/centromere-prediction` (`git fetch` + SHA diff, no API call). Verdicts: `UP-TO-DATE` / `LOCAL-AHEAD` (exit 0); `BEHIND-BY-N` (1: rsync from the source checkout, then `diff -rq`); `OFFLINE` / `NO-ORIGIN` (2, informational). Implementation: `bin/skill-update-check.py`.

## Glossary

- **Monomer**: consensus tandem repeat unit (default 50–250 bp, up to 2000 with `--period_max`) discovered from reads and clustered with cd-hit-est.
- **HOR**: Higher-Order Repeat — hierarchical monomer arrangement typical of vertebrate centromeres; visible as periodic dominant-monomer structure in the decomposition.
- **Density window**: 100 kb interval scored by fraction of bases covered by candidate repeats (column 7 of `centromere_density_100kb.bed`).
- **Candidate region**: the densest `region_window` (default 1 Mb) per scaffold (`candidate_regions.tsv`) — a *candidate*, not a proven centromere. TRF reports a clean array as ONE long repeat, so density 1.0 with 1–3 repeats is a real array, not an artifact.
- **`publishDir overwrite: false`**: stale outputs from earlier runs are silently kept; always move old `results/` aside before re-running.

## Verification

- [ ] Preflight performed; verdict ≥ `GO-WITH-WARNINGS`.
- [ ] Run completed; `results/centromere/centromere_summary.tsv` present.
- [ ] QC performed; `centromere_report.md` exists with verdict; a negative verdict was re-checked with a wide period range.

## Handoff pointers

After routing the user to a phase, state explicitly which sub-skill was invoked and what file it produces.

## Provenance

- Rebuilt with the [skill-forge](https://github.com/cheahhl814/skill-forge) meta-skill (`nf-centromere-agents.zip` as source; BettaMt pattern).
- Worked examples: PR001 *Parosphromenus rubrimontis* (v1.0.0: 500 Mb, 173 contigs, < 1 h on 24 cpus / 96 GB) and *Betta coccina* BC001 (v1.1.0: 478 Mb, 33 scaffolds; a 685 bp satellite array found only with `--period_max 2000`, on 21 of 23 scaffolds ≥ 5 Mb).
