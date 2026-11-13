# centromere-prediction

Centromere annotation & analysis for ONT-based fish genome assemblies (Betta_v3 output).

Discovers tandem-repeat monomers from sequencing reads (TideHunter + cd-hit-est), maps
them onto the assembly (TRF + bedtools density), picks the highest-density region per
scaffold, decomposes each region into monomers (stringDecomposer), and produces a
per-scaffold summary — packaged as a BettaMt-convention agent skill wrapping a thin
Nextflow DSL2 runner.

> [!NOTE]
> This skill is a rebuild of the standalone `centromere-prediction` pipeline (5 latent
> Betta_v3-era bugs fixed, see the orchestrator `SKILL.md` §4). **v1.1.0** (first real
> *Betta coccina* run) makes the period range a parameter, bounds `stringdecomposer`
> memory, publishes the real decomposition table (v1.0.0 published only its log), and adds
> satellite mapping, a telomere cross-check and evidence figures.

## Contents

- [Installation](#installation)
- [Usage](#usage)
- [Pipeline overview](#pipeline-overview)
- [Tools](#tools)
- [Update check](#update-check)
- [Repository layout](#repository-layout)
- [Hard guarantees](#hard-guarantees)
- [Provenance](#provenance)
- [License](#license)

## Installation

**For an AI agent**, prompt it with a natural-language trigger:

- *"find centromeres on my fish assembly"* / *"check centromere results"*
- *"run centromere analysis with TideHunter and stringdecomposer"*
- the sub-skill names: `nf-centromere-run-preflight`, `nf-centromere-pipeline-run`, `nf-centromere-results-qc`, `nf-centromere-run-debug`

**Manually** (self-contained pixi manifest; tool env auto-installs on first use):

```bash
git clone https://github.com/cheahhl814/centromere-prediction.git
cd centromere-prediction
pixi lock --dry-run           # verify the 7-tool env solves
pixi run -e nextflow nextflow -version   # optional isolated runner env
```

Nextflow ≥ 23.10 required — on PATH, or materialize the isolated env:
`pixi run -e nextflow nextflow`.

**Agent-skill wiring** (hard deps, declared in `SKILL.md` frontmatter `skills:`):
upstream `read-qc-trimming` (trimmed nuclear reads input) and `chartz`
(post-QC visualisation of centromere evidence — e.g. `n_repeats` per scaffold from
`centromere_summary.tsv`, or repeat-density profiles). No database download:
all 6 tools consume only user-supplied files; the only network activity is pixi
package installation on first use (tool binaries, not reference databases).

## Usage

### Natural-language triggers

> "Find candidate centromeres on this assembly" ·
> "Run the centromere pipeline on my trimmed reads" ·
> "Check my centromere results" ·
> "The centromere run failed — diagnose it" ·
> "Loosen the monomer period range"

### Manual — the four phases

The agent sub-skills are the source of truth; these commands map 1:1:

| Phase | Command (run from your analysis run dir) |
|---|---|
| Preflight | see `preflight/run-preflight/SKILL.md` (input + env + disk validation → `params.json`) |
| Run | `pixi run -e nextflow nextflow run <skill-root>/runners/nextflow-runner/main.nf -profile local -params-file centromere_params.json` |
| Preview (dry-run) | same command + `-preview` (Nextflow ≥ 25.x; older versions use `-dry-run`) |
| QC + report | see `qc/results-qc/SKILL.md` → writes `centromere_report.md` |
| Debug | see `debug/run-debug/SKILL.md` |

Environment task stubs exist as reminders: `pixi run preflight|run|qc|debug`,
plus `pixi run lock-check` (dry-run solve of every env) and
`pixi run update-check`.

A worked-example params file ships at
`runners/nextflow-runner/params.example.json` — replace its paths with your own
assembly/reads before running.

## Pipeline overview

| Process | Output | Purpose |
|---|---|---|
| `tidehunter_cluster` | `tidehunter/consensus_clustered.fasta`, `monomers_top.fasta`, `monomer_ranking.tsv` | tandem-repeat monomers from reads, clustered at 80%, ranked by cluster size, top `--max_monomers` kept |
| `trf_density` | `trf/all_repeats.bed`, `trf/centromere_candidate_repeats.bed`, `trf/centromere_density_100kb.bed`, `trf/genome.sizes` | repeats in `--period_min`–`--period_max` bp with ≥ `--min_copies`, density per 100 kb window |
| `identify_regions` | `regions/candidate_regions.tsv` | densest `--region_window` (default 1 Mb) per scaffold |
| `stringdecomposer` | `stringdecomposer/<scaffold>/decomposition.tsv` (+ log, region fasta) | monomer decomposition of each candidate region against the capped monomers |
| `centromere_summary` | `centromere_summary.tsv` | per-scaffold repeat span + count of TRF entries |

Parameters (all optional): `period_min` 50, `period_max` 250 (TRF caps at 2000), `min_copies` 5,
`region_window` 1000000 (0 = legacy merged top-10 span), `max_monomers` 100 (0 = all).
Run a **second, wide pass** (`period_max: 2000`) before accepting a negative result: TRF stops
at `period_max` itself, so a post-filter cannot recover longer satellites.

QC helpers (run after the pipeline; see `qc/results-qc/SKILL.md` steps 3 and 6–8), all under `bin/`:
`decomp_stats.py` (per-scaffold dominant-monomer statistics), `satellite_map.py` (BLASTn the dominant
monomer genome-wide and merge hits into arrays), `telomere_crosscheck.py` (tidk; internal vs capped vs
uncapped near-end arrays), `evidence_tracks.py` (per-scaffold density/satellite/telomere figures via chartz).

The interpretation rules (healthy ranges, red flags, verdicts) live in
`qc/results-qc/SKILL.md` — the single source of truth for reading results.

## Tools

All 8 bioinformatics tools co-solve in the default env (verified by `pixi search`
+ joint solve at build time); `nextflow` lives in the isolated `nextflow` env
(`no-default-feature`), materialized only on demand.

| Tool | Version pin | Role |
|---|---|---|
| TideHunter | ≥1.5.5 | monomer discovery from reads |
| TRF | ≥4.09.1 | tandem repeat mapping on assembly |
| cd-hit | ≥4.8.1,<5 | monomer consensus clustering |
| samtools | ≥1.19 | faidx extraction + genome sizes |
| bedtools | ≥2.31 | density windows + coverage + summary |
| stringdecomposer | ≥1.1.2 | monomer decomposition of candidate regions |
| BLAST+ | ≥2.14 | map the dominant satellite genome-wide (`satellite_map.py`) |
| tidk | ≥0.2.6 | telomere-repeat windows (`telomere_crosscheck.py`) |
| nextflow | ≥23.10 (isolated env) | pipeline runtime (DSL2) |

Offline flag/usage references live in `docs-corpus/<tool>/` — the agent reads
those instead of guessing CLI flags.

## Update check

```bash
pixi run update-check
# UP-TO-DATE (0) | LOCAL-AHEAD (0) | BEHIND-BY-N (1) → rsync from source checkout |
# OFFLINE (2) | NO-ORIGIN (2) — informational
```

Implemented in `bin/skill-update-check.py` (stdlib-only, git-SHA diff via `git fetch`,
no GitHub API).

## Repository layout

```
centromere-prediction/
├── SKILL.md                    # orchestrator router (start reading here for agents)
├── README.md                   # this file
├── pixi.toml                   # one manifest: 8 tools + isolated nextflow env
├── params.json                 # skill-forge render input (tracked)
├── bin/
│   ├── extract_candidate_region.sh   # densest window per scaffold (window 0 = legacy merged top-10)
│   ├── rank_monomers.py              # cd-hit clusters → top-N by size + ranking table
│   ├── decomp_stats.py               # per-scaffold dominant-monomer stats
│   ├── satellite_map.py              # BLASTn a monomer genome-wide → arrays
│   ├── telomere_crosscheck.py        # tidk: internal vs capped vs uncapped arrays
│   ├── evidence_tracks.py            # chartz per-scaffold evidence figures
│   └── skill-update-check.py         # self-update check
├── preflight/run-preflight/    # Phase 1: input/env validation
├── run/pipeline-run/           # Phase 2: run the workflow
├── qc/results-qc/              # Phase 3: interpret → centromere_report.md
├── debug/run-debug/            # failure interpreter (signature library)
├── docs-corpus/<tool>/         # offline tool docs (help captures / repo READMEs)
├── test-data/                  # deterministic fixtures + MANIFEST.tsv (+ satellite fixture, MANIFEST_sat.tsv)
└── runners/nextflow-runner/    # the thin DSL2 runner (main.nf, modules, config)
```

## Hard guarantees

- The pipeline ran end-to-end on two real assemblies: PR001 (v1.0.0, 173 scaffolds) and
  *B. coccina* BC001 (v1.1.0, 33 scaffolds; the stock v1.0.0 step OOM-killed on it, which
  motivated the v1.1.0 changes). `-profile test -stub-run` (empty-signal fixtures) and
  `-profile test_sat` (a real run on a satellite fixture) are the regression canaries.
- The default period range (50–250 bp, ≥5 copies) is vertebrate-tuned. It is a parameter, and
  a negative result must be re-checked at `period_max: 2000` before it is believed.
- `publishDir` uses `overwrite: false` — stale-run contamination is guarded by the
  preflight stale-results check and surfaced in run/debug instructions.
- Deterministic test fixtures (`test-data/`, fixed seed, sha256 manifest) back the
  structural smoke tests.

## Provenance

Rebuilt from `centromere-prediction.zip` with the
[skill-forge](https://github.com/cheahhl814/skill-forge) meta-skill (BettaMt
pattern, cf. [BettaMt-agents](https://github.com/cheahhl814/BettaMt-agents),
[Betta-WGS-agent](https://github.com/cheahhl814/Betta-WGS-agent)). Worked example:
*Parosphromenus rubrimontis* (PR001) — TideHunter ~13 min + TRF ~26 min in parallel,
full pipeline < 1 h on 24 cpus / 96 GB.

## License

See repository settings / LICENSE when published.