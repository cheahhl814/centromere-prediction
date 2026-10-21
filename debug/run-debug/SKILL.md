---
name: nf-centromere-run-debug
description: 'Failure interpreter for nf-centromere-agents: match a failed centromere run (stderr, run log, or weird outputs) against a signature library of known failure patterns — channel semantics, publishDir overwrite:false stale files, -resume gotchas, empty TRF/TideHunter outputs, tool-env gaps, and legacy Betta_v3 quirks — and recommend a concrete fix. Use when a centromere run fails, produces empty or partial outputs, or behaves suspiciously (e.g. stringdecomposer ran for only one scaffold).'
version: 1.0.0
updated: "2026-10-21"
triggers:
  - "centromere run failed"
  - "centromere debug"
  - "diagnose centromere"
---

# Sub-skill: run-debug (Failure Interpreter)

## 0. Inputs / Outputs contract

| Direction | Artifact | Notes |
| --- | --- | --- |
| Inputs | failing stderr / `centromere_run.log` tail / nextflow `.trace` or work-dir dir for the failed process | |
| Outputs | `<RUN_DIR>/debug-report.md` with diagnosis + fix (verdict `RECOVERABLE` or `BLOCKED`) | |

## Procedure

1. Read the actual error — `tail -50 "$RUN_DIR/centromere_run.log"`, or `cd "$RUN_DIR" && tree .nextflow/log/` and open the newest `.command.err` under the failed process's work dir. Never conclude from memory.
2. Match against the signature library below. First match wins.
3. If no match: short investigation loop — re-read the matching `docs-corpus/<tool>/` entry for exact flags; check tool versions (`pixi run --manifest-path <root>/pixi.toml <tool> --version`); re-run the single process manually with `nextflow run ... -process X` or the bash recipe.
4. Write `debug-report.md`: signature matched, root cause, concrete fix (command + tool/version), verdict.
5. If the fix required a pipeline change, record it: patch `runners/nextflow-runner/modules/centromere.nf` and bump the skill version in the frontmatter — then note in `debug-report.md` that the fix is in-tree (never fix out-of-tree per-run copies).

## Signature library

| Signature in stderr / log | Likely cause | Suggested fix |
| --- | --- | --- |
| `centromere` process fails; `.command.err` shows `command not found` for any of TideHunter/TRF/cd-hit/samtools/bedtools/stringdecomposer | pixi env not materialized for the runner's manifest path | `pixi run --manifest-path <skill-root>/pixi.toml <tool> --version`; first materialization needs network; then `nextflow run ... -resume`. |
| Only `stringdecomposer/<one-scaffold>/` produced output while `candidate_regions.tsv` has many rows | pre-`143a5dc` queue-channel bug (assembly/monomers passed as queue channels) | Verify runner has `Channel.value(file(...))` in `main.nf`; trash stale `results/centromere/stringdecomposer/` (publishDir `overwrite: false` keeps old files silently); run fresh. |
| `.command.err`: `Error: Could not recognize the command: preview` (old NF versions inverse) or `-dry-run` unknown | Nextflow version mismatch on the preview flag | Nextflow ≥25 uses `-preview`; older uses `-dry-run`. Check `nextflow -version`. |
| `consensus_clustered.fasta` missing from `results/centromere/tidehunter/` but present in work dirs | legacy publish-pattern typo (`consensus_clustered.fa`) — this repo fixed it | Verify module pattern is `.fasta`; for old runs, salvage: `find .nextflow/work -name consensus_clustered.fasta`. |
| `candidate_regions.tsv` missing from `results/centromere/regions/` | legacy: `identify_regions` had no publishDir | This runner publishes it; if missing in an old run, copy from the work dir. |
| `.command.err`: `bedtools: command not found` inside `centromere_summary` | legacy: bare `bedtools` outside the env | This runner wraps it in `pixi run`; for old runs install bedtools globally or wrap manually. |
| `centromere_density_100kb.bed` empty | TRF produced no output; the module's `\|\| true` swallowed the error | Re-run TRF standalone: `trf <assembly> 2 7 7 80 10 50 250 -f -d -m -ngs -h`; inspect its stderr directly (without `\|\| true`). |
| `consensus_clustered.fasta` empty (file exists, 0 bytes) | no tandem-repeat signal in reads | Verify reads are long (>10 kb ideal); check `-p 50 -P 250`; suspect mtDNA-polluted or heavily-shortened reads. |
| TideHunter errors on `.gz` input | tool doesn't accept compressed input | Decompress first (`unpigz -c in.fastq.gz > in.fastq`); keep decompressed file outside the run-dir git scope. |
| `tidehunter_cluster` crashes with a process arity error (2 inputs declared, 1 supplied) | legacy: unused `val clust_id` input | This runner removed it; do not resurrect it. |
| `-resume` replays everything but output didn't change | resumed the wrong session — `.nextflow/history` last entry was a `-preview`/aborted run | Run fresh: `nextflow run ... -params-file ...` (no `-resume`) after moving old results aside. |
| `centromere_summary.tsv` missing/empty, everything upstream published | no repeats passed the period/copies filter | QC sub-skill: loosen period range (post-process) → verdict `monomer_filter_too_strict`; check `all_repeats.bed` first. |
| `Permission denied` in work dir | work dir on a read-only mount or wrong owner | Point `-work-dir` to a writable path on the same filesystem as `$RUN_DIR`. |
| Sudden OOM / SIGKILL on `trf_density` (huge genome) | 100 kb windows over a very large genome with dense repeats; 96 GB ceiling | Reduce scope: pre-filter scaffolds, or request more memory via `nextflow.config` `withLabel: 'centromere'`. |

## Legacy Betta_v3 results note

Old results were produced by the buggy in-Betta_v3 module. For interpretation of legacy dirs: treat a single-scaffold `stringdecomposer/` as **incomplete**, not as "only one scaffold has a centromere" — re-run the fixed runner before drawing biological conclusions.

## Handoff

After writing `debug-report.md`:

- **`RECOVERABLE`** → back to `run/pipeline-run` (usually fresh-run after moving stale outputs).
- **`BLOCKED`** → hand to the user with the unresolvable evidence (e.g. corrupted inputs, missing upstream data).