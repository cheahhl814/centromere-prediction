---
name: centromere-prediction-pipeline-run
description: 'Phase 2 of centromere-prediction: run the 5-process centromere Nextflow workflow (TideHunter monomer discovery, TRF density, identify_regions, stringdecomposer, centromere_summary) with a dry-run preview, explicit user confirmation, background-execution pattern, and the -resume / publishDir gotchas. Use when preflight passed and the centromere run should be started, resumed, or monitored.'
version: 1.0.0
updated: "2026-10-21"
triggers:
  - "run centromere pipeline"
  - "start centromere run"
  - "resume centromere"
---

# Sub-skill: pipeline-run (Phase 2 — The Workflow)

## 0. Inputs / Outputs contract

| Direction | Artifact | Notes |
| --- | --- | --- |
| Inputs | `<RUN_DIR>/params.json` + `preflight_passed.md` (verdict ≥ GO-WITH-WARNINGS) | from run-preflight |
| Outputs | `<RUN_DIR>/results/centromere/**` (6 artifact kinds) | consumed by results-qc |
| Outputs | `<RUN_DIR>/centromere_run.log` (`results/pipeline_info/` trace/report/timeline) | |

Refuse to run if `preflight_passed.md` is missing or its verdict is `NO-GO` — tell the user which preflight check to fix.

## Pipeline map (read this to understand what runs)

| Process | Input | Output | Does |
|---|---|---|---|
| `tidehunter_cluster` | reads (fastq) | `tidehunter/consensus_clustered.fasta` | Finds tandem-repeat monomers (period 50–250, `-c 5`), clusters at 80% identity (`cd-hit-est -c 0.8`) |
| `trf_density` | assembly fasta | `trf/centromere_candidate_repeats.bed`, `trf/centromere_density_100kb.bed`, `trf/genome.sizes` | Runs TRF, keeps repeats with period 50–250 bp and ≥5 copies, computes density in 100 kb windows |
| `identify_regions` | density bed | `regions/candidate_regions.tsv` | Top-10 density windows per scaffold merged into one candidate region (`bin/extract_candidate_region.sh`) |
| `stringdecomposer` | tuple(scaffold, region, monomers) + assembly **value** channel | `stringdecomposer/<scaffold>/decomposition.tsv`, `<scaffold>_centromere.fa` | Decomposes each candidate region into monomers |
| `centromere_summary` | candidate repeats bed | `centromere_summary.tsv` | Per-scaffold start/end/count of candidate repeats |

`tidehunter_cluster` and `trf_density` run in parallel; the rest is sequential. No GPU needed.

## Procedure

### 1. Build the params file (from the preflight contract)

```bash
cd "$RUN_DIR"
python3 - <<'EOF'
import json
p = json.load(open("params.json"))
json.dump({"reads": p["reads"], "assembly": p["assembly"]}, open("centromere_params.json", "w"), indent=2)
EOF
```

### 2. Dry-run preview first: `NXF` 25.x uses `-preview`, not `-dry-run`

```bash
NF="nextflow"
pixi run -e nextflow nextflow -version >/dev/null 2>&1 && NF="pixi exec -e nextflow nextflow"
$NF run "$SKILL_ROOT/runners/nextflow-runner/main.nf" -profile local -params-file centromere_params.json -preview
```

### 3. Ask: "Start the run? (yes / no)" — never run without explicit confirmation

### 4. Launch (background, interruption-safe)

```bash
cd "$RUN_DIR"
TERM=xterm setsid nohup $NF run "$SKILL_ROOT/runners/nextflow-runner/main.nf" \
    -profile local -params-file centromere_params.json \
    > centromere_run.log 2>&1 < /dev/null &
echo "launched; log: $RUN_DIR/centromere_run.log"
```

Progress: `tail -30 "$RUN_DIR/centromere_run.log"`.

**Observed timings (PR001, 500 Mb assembly, 173 contigs, 24 cpus / 96 GB):** TideHunter ~13 min and TRF ~26 min (parallel), identify/summary instant, stringdecomposer = 173 short tasks (seconds each). Full pipeline **under an hour** — when in doubt, run fresh rather than fighting `-resume`.

### 5. Resume after failure

```bash
$NF run "$SKILL_ROOT/runners/nextflow-runner/main.nf" -profile local -params-file centromere_params.json -resume
```

**Resume gotchas:**

- `-resume` takes **no argument** — it resumes from the **last session in `.nextflow/history`**. If the last session was a `-preview` dry-run or an aborted/empty run, resume starts from scratch.
- The modules publish with `overwrite: false` — a resumed/fresh run **silently keeps stale files** from earlier runs. Before re-running, move the old `results/centromere/` to trash.

### 6. Manual single-stage re-run (bash recipe — debugging only)

The bash recipes mirror the modules and share the skill-root pixi manifest. Example (TRF standalone, as used by the debug sub-skill):

```bash
pixi run --manifest-path "$SKILL_ROOT/pixi.toml" \
    trf <assembly> 2 7 7 80 10 50 250 -f -d -m -ngs -h > assembly.trf.dat 2>/dev/null
```

TideHunter does **not** accept compressed input — decompress first (`unpigz -c in.fastq.gz > in.fastq`).

### Integration

This pipeline is intentionally standalone (the module was removed from Betta_v3). Run it after each species' WGS assembly, using that run's trimmed reads + polished assembly.

## Signature library

| Signature in run log | Likely cause | Fix |
| --- | --- | --- |
| `Missing required parameters: provide --reads and --assembly` | params file not passed or keys wrong | Re-generate `centromere_params.json` from `params.json`; verify keys. |
| `Channel` error about single-value queue channel at `stringdecomposer` | running the pre-`143a5dc` Betta_v3 module | Use this repo's runner (value-channel fix); verify via `grep -n "Channel.value" runners/nextflow-runner/main.nf`. |
| Workflow stops within seconds, log ends with `--preview` output only | you ran `-preview` instead of the real run | Run again without `-preview` after confirmation. |
| All steps `CACHED` but no new outputs | stale results + `overwrite: false` | Move `results/centromere/` aside; run fresh (`-resume` will replay cached work). |
| `Unknown option: -dry-run` | Nextflow 25.x renamed the flag | Use `-preview`. |
| `Process ... failed: Bedtools` etc. inside a process | tool env not materialized for that process | `pixi run --manifest-path $SKILL_ROOT/pixi.toml <tool> --version` standalone; re-run with `-resume`. |
| `<tool>: command not found` inside a task where the same command works manually (and nextflow was launched via `pixi run -e nextflow`) | nested-pixi env resolution: `PIXI_*` activation vars inherited from the parent pixi run make the task's `pixi run --manifest-path` silently resolve the active (nextflow) env | This repo's modules already prefix every pixi call with `env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE`; verify with `grep -c "env -u PIXI" runners/nextflow-runner/modules/centromere.nf` or launch nextflow from a non-pixi (PATH) install |
| Script compile error `Statements cannot be mixed with script declarations` at `main.nf` top-level `if` | Nextflow ≥26 hard-errors on bare statements at script top level | Move the params validation inside the `workflow` block (this repo's `main.nf` already does) |
| `No such property: <input var> for class: ProcessDslV1` at a `publishDir` line | strict-DSL parser v2 evaluates bare-string publishDir at process-definition time, before inputs exist | Wrap the path in a closure: `publishDir path: { ".../${scaffold}" }, mode: ..., ...` (this repo's `stringdecomposer` already does) |