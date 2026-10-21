---
name: nf-centromere-run-preflight
description: 'Phase 1 of nf-centromere-agents: validate centromere run inputs (reads FASTQ/FASTA + assembly FASTA), verify the pixi tool env solves, check the helper script and Nextflow availability, and write the run params.json + preflight verdict (GO / GO-WITH-WARNINGS / NO-GO). Use when starting a centromere analysis, after locating trimmed nuclear reads and a polished assembly, or before invoking the centromere run sub-skill.'
version: 1.0.0
updated: "2026-10-21"
triggers:
  - "centromere preflight"
  - "validate centromere inputs"
  - "prepare centromere run"
---

# Sub-skill: run-preflight (Phase 1 — The Audit)

> Validates the centromere run's inputs and environment. Consumes nothing but user-provided paths; produces `params.json` + `preflight_passed.md` in the **analysis run dir**.

## 0. Inputs / Outputs contract

| Direction | Artifact | Notes |
| --- | --- | --- |
| Inputs | reads path (fastq/fasta, `.gz` ok for reads), assembly path (fasta) | from the user or the upstream WGS run |
| Inputs | `$SKILL_ROOT/runners/nextflow-runner/` | the pipeline to be run |
| Outputs | `<RUN_DIR>/params.json` (machine contract for the run) | |
| Outputs | `<RUN_DIR>/preflight_passed.md` (verdict gate for run) | |

Upstream skills that produce these inputs (declared as `skills:` deps in the master frontmatter):

- **`read-qc-trimming`** → trimmed nuclear reads (e.g. `qc/nf-qc-run/results/trim/<species>_trimmed.fastq.gz`). If the user offers raw reads without trimming, route there first; mtDNA/organellar tandem repeats in untrimmed raw reads pollute the monomer set. (The assembly input comes from whatever upstream WGS/pipeline run the user has — no assembly-building skill is hard-wired.)
- **`chartz`** → after the QC sub-skill's interpretation, when the user wants visualisation of the centromere evidence (e.g. a bar chart of `n_repeats` per scaffold from `centromere_summary.tsv`, or repeat-density profile plots from `centromere_density_100kb.bed`). Hard-wired: QC always offers it; every plotted value is a named column of the existing results files.

Downstream skills (prose only — not hard wiring):
- `ena-submit` — if the user wants to publish the assembly, not the centromere outputs.

Refuse to proceed (and route the user back) if the analysis `RUN_DIR` is not writable or both input paths are unresolved.

## Ask-user stop points

- **SP1 — Input provenance ambiguous**: reads could be raw ONT reads vs trimmed nuclear reads vs mtDNA-unfiltered. Evidence + Recommend + Options: recommend trimmed nuclear reads (`qc/nf-qc-run/results/trim/<species>_trimmed.fastq.gz` style); accept raw reads only with a warning that mtDNA/organellar tandem repeats pollute the monomer set.
- **SP2 — Assembly vintage**: polished+scaffolded assembly preferred; a draft/contig-level assembly fragments candidates (see QC red flags). Warn, allow proceed.

## Procedure

```bash
RUN_DIR="${RUN_DIR:-$PWD}"          # analysis run dir — where results land
SKILL_ROOT=<path to this skill>     # e.g. ~/.agents/skills/nf-centromere-agents

# 1. Inputs exist and are non-empty
for f in "<reads.fastq.gz>" "<assembly.fasta>"; do
    test -s "$f" || { echo "PREFLIGHT FAIL: missing/empty input: $f"; exit 1; }
done

# 2. Read shape (sanity): gzip integrity + record count for reads
zcat "<reads.fastq.gz>" | head -4            # well-formed FASTQ? (fastq) — or grep '>' for fasta
# Assembly: >1 record, plausible N50 context (run seqkit stats if available — optional)

# 3. Tool env solves
pixi lock --dry-run --manifest-path "$SKILL_ROOT/pixi.toml" || {
    echo "PREFLIGHT FAIL: pixi env unsolvable"; exit 1; }
# nextflow runtime: on PATH, or materialize the isolated env
command -v nextflow >/dev/null 2>&1 || pixi run -e nextflow nextflow -version

# 4. Helper script + pipeline entry present
test -x "$SKILL_ROOT/bin/extract_candidate_region.sh"
test -f "$SKILL_ROOT/runners/nextflow-runner/main.nf"

# 5. Disk: TRF intermediates can reach a few GB on large genomes
FREE_GB=$(df --output=avail -BG "$RUN_DIR" | tail -1 | tr -dc '0-9')
test "$FREE_GB" -ge 5 || echo "WARN: only ${FREE_GB} GB free"

# 6. Stale-results guard: publishDir uses overwrite:false
if [ -d "$RUN_DIR/results/centromere" ]; then
    echo "WARN: stale results/centromere/ exists — move it aside before running or old files will be kept silently."
fi
```

Then write `params.json` (run contract consumed by the run sub-skill):

```json
{
  "reads": "<absolute path to reads>",
  "assembly": "<absolute path to assembly>",
  "profile": "local",
  "skill_root": "<SKILL_ROOT>"
}
```

...and `preflight_passed.md` with a verdict table (inputs / env / runtime / disk / stale-results → `GO` | `GO-WITH-WARNINGS` | `NO-GO`). Never run with a `NO-GO`.

## Signature library

| Signature in preflight | Likely cause | Fix |
| --- | --- | --- |
| `pixi lock --dry-run` fails with resolution error | conda-forge/bioconda index out of date or offline | Retry once; run `pixi --version` and check network; report verbatim error. |
| `nextflow: command not found` after env check | nextflow not on PATH and isolated env not materialized | `pixi run -e nextflow nextflow -version` (first use downloads java + nextflow). |
| Reads are FASTA but named `.fastq.gz` | mislabeled input | Inspect first bytes (`zcat \| head -1`); fix the params path, not the file. |
| Assembly has 0 or 1 record | wrong file (contigs vs scaffolds confusion) | Re-locate upstream `ntLink_scaffolds_polish.fasta`; ask user (SP2). |
| `zcat` fails mid-file (unexpected EOF) | truncated download | Re-fetch the reads from upstream QC outputs; do not proceed. |