# centromere-prediction — Centromere Annotation & Analysis

[![Version](https://img.shields.io/badge/version-1.1.0-blue)](#-installation)
[![Type](https://img.shields.io/badge/type-agent%20skill-blueviolet)](#-installation)
[![Built with](https://img.shields.io/badge/built%20with-bioinfo--skill--creator-orange)](https://github.com/cheahhl814/skill-forge)

[![Shape](https://img.shields.io/badge/shape-multi--tool%20pipeline-8A2BE2)](#-pipeline-overview)

Centromere annotation and analysis for ONT-based fish genome assemblies. Discovers tandem-repeat monomers from reads (TideHunter + cd-hit-est), maps them onto the assembly (TRF + bedtools density windows), picks the densest region per scaffold, decomposes each region into monomers (stringDecomposer), and produces a per-scaffold summary. A QC phase then maps the dominant satellite genome-wide (BLASTn), cross-checks it against telomere repeats (tidk), plots the evidence, and writes a verdict report. Packaged as a BettaMt-convention agent skill wrapping a thin Nextflow DSL2 runner. Outputs are **candidate** centromeres, not proven ones.

**Repository**: https://github.com/cheahhl814/centromere-prediction

> [!NOTE]
> Current version: **v1.1.0** (updated 2026-11-03). v1.1.0 makes the monomer period range a parameter, bounds stringDecomposer memory, publishes the real decomposition table (v1.0.0 published only its log), and adds satellite mapping, a telomere cross-check and per-scaffold evidence figures.

## Contents

- [Installation](#-installation)
- [Usage](#-usage)
- [Pipeline overview](#-pipeline-overview)
- [Tools](#-tools)
- [Update check](#-update-check)
- [Depends on](#-depends-on)
- [Repository layout](#-repository-layout)
- [Hard guarantees](#-hard-guarantees)
- [Nextflow runner](#-nextflow-runner)
- [Provenance](#-provenance)

## 🚀 Installation

This is an **agent skill**, not a user-facing library. The recommended install path is to let your AI agent import it.

**Option A — give your agent this prompt (recommended):**

```text
Install the centromere-prediction skill from
https://github.com/cheahhl814/centromere-prediction —
clone it into your agent's skills directory (the path your agent watches
for skills) and run `pixi install` from that directory. Then read the
skill's SKILL.md to understand its phases and the Input/Output contract
for each. Confirm when the environment is ready.
```

**Option B — manual install:**

```bash
git clone https://github.com/cheahhl814/centromere-prediction.git
cd centromere-prediction
pixi install
```

> [!TIP]
> `pixi install` resolves every pinned tool from `pixi.toml` (channels: conda-forge, bioconda) into an isolated `.pixi/` environment — no system-wide installs, no version conflicts with other skills. Nextflow lives in a separate on-demand env: `pixi run -e nextflow nextflow -version`.

## 💡 Usage

The skill is designed to be driven by an AI agent: the agent reads the master `SKILL.md`, detects the current phase from the filesystem state of your run directory, and routes to the right sub-skill.

### Natural-language prompts that trigger the skill

```text
find centromeres on my fish assembly
run centromere analysis with TideHunter and stringdecomposer
check centromere results
the centromere run failed — diagnose it
widen the monomer period range
```

### Manual phase execution

If you prefer to drive the phases yourself (or want to re-run a single phase):

```bash
# Set your run directory (all phase artifacts land here)
RUN_DIR=/path/to/run-dir

pixi run preflight   # Phase 1 — validate inputs, write params.json (see preflight/run-preflight/SKILL.md)
pixi run run         # Phase 2 — execute the workflow (see run/pipeline-run/SKILL.md)
pixi run qc          # Phase 3 — interpret results, write centromere_report.md (see qc/results-qc/SKILL.md)
pixi run debug       # On failure — interpret stderr via the signature library
pixi run lock-check  # Dry-run solve of every environment
```

The `preflight`, `run`, `qc` and `debug` tasks print a pointer to the matching sub-skill; the real commands live there. The pipeline itself is launched with:

```bash
pixi run -e nextflow nextflow run <skill-root>/runners/nextflow-runner/main.nf \
  -profile local -params-file centromere_params.json
```

A worked-example params file ships at `runners/nextflow-runner/params.example.json` — replace its paths with your own assembly and reads.

> [!IMPORTANT]
> Each phase has an explicit **Inputs/Outputs contract** at the top of its sub-skill `SKILL.md`. If the upstream artifact is missing (e.g. you run `run` before `preflight` passed), the sub-skill refuses to proceed and tells you which phase to run first. Phases are gated on purpose — don't skip them.
>
> **For a new species, run twice:** the default 50–250 bp range first, then a wide pass with `period_max: 2000` in a separate run dir. A negative result from the default range is provisional, because TRF stops at `period_max` itself and a post-filter cannot recover longer satellites.

## 🔬 Pipeline overview

The analysis follows a phased evidence chain. Each phase consumes the artifacts of the previous phase and gates progression with a Go/No-Go check.

| # | Phase | Goal | Sub-skill | Artifact produced |
|:--|:------|:-----|:----------|:------------------|
| **1** | **Preflight** | Validate reads, assembly, tool env and disk; pick the read subsample and memory-safe parameters | `preflight/run-preflight/` | `$RUN_DIR/params.json`, `$RUN_DIR/preflight_passed.md` |
| **2** | **Run** | Execute the 5-process Nextflow workflow | `run/pipeline-run/` | `results/centromere/centromere_summary.tsv` |
| **3** | **QC** | Interpret outputs, map and cross-check the satellite, plot evidence, give a verdict | `qc/results-qc/` | `$RUN_DIR/centromere_report.md` |
| — | **Debug** | Match a failure to the signature library | `debug/run-debug/` | fix recommendation |

The Nextflow workflow (phase 2) has five processes:

| Process | Output | Purpose |
|:--------|:-------|:--------|
| `tidehunter_cluster` | `tidehunter/consensus_clustered.fasta`, `monomers_top.fasta`, `monomer_ranking.tsv` | Tandem-repeat monomers from reads, clustered at 80%, ranked by cluster size, top `--max_monomers` kept |
| `trf_density` | `trf/all_repeats.bed`, `trf/centromere_candidate_repeats.bed`, `trf/centromere_density_100kb.bed`, `trf/genome.sizes` | Repeats in `--period_min`–`--period_max` bp with ≥ `--min_copies`, density per 100 kb window |
| `identify_regions` | `regions/candidate_regions.tsv` | Densest `--region_window` (default 1 Mb) per scaffold |
| `stringdecomposer` | `stringdecomposer/<scaffold>/decomposition.tsv` (+ log, region fasta) | Monomer decomposition of each candidate region against the capped monomers |
| `centromere_summary` | `centromere_summary.tsv` | Per-scaffold repeat span and count of TRF entries |

All parameters are optional:

| Parameter | Default | Meaning |
|:----------|:--------|:--------|
| `period_min` / `period_max` | 50 / 250 | Monomer period range in bp (TRF caps at 2000) |
| `min_copies` | 5 | Minimum copies for a TRF candidate |
| `region_window` | 1000000 | Size of the densest window to decompose per scaffold (0 = legacy merged top-10 span) |
| `max_monomers` | 100 | Most abundant cd-hit clusters passed to stringDecomposer (0 = all) |

QC helpers, run after the pipeline (see `qc/results-qc/SKILL.md`), all under `bin/`:

| Helper | Purpose |
|:-------|:--------|
| `decomp_stats.py` | Per-scaffold dominant-monomer statistics |
| `satellite_map.py` | BLASTn the dominant monomer genome-wide and merge hits into arrays |
| `telomere_crosscheck.py` | tidk telomere repeats: internal vs capped vs uncapped near-end arrays |
| `evidence_tracks.py` | Per-scaffold density / satellite / telomere figures via chartz |

> [!TIP]
> Read each sub-skill's `SKILL.md` for the full procedure, its ask-user stop points, and its signature library (stderr pattern → cause → fix). The interpretation rules (healthy ranges, red flags, verdicts) live in `qc/results-qc/SKILL.md`.

## 🧰 Tools

All tools are resolved from conda-forge/bioconda via the pinned `pixi.toml`.

| Tool | Version | Role |
|:-----|:--------|:-----|
| `tidehunter` | ≥1.5.5 | Monomer discovery from reads |
| `trf` | ≥4.09.1 | Tandem repeat mapping on the assembly |
| `cd-hit` | ≥4.8.1,<5 | Monomer consensus clustering |
| `samtools` | ≥1.19 | faidx extraction + genome sizes |
| `bedtools` | ≥2.31 | Density windows, coverage and summary |
| `stringdecomposer` | ≥1.1.2 | Monomer decomposition of candidate regions |
| `blast` | ≥2.14 | Map the dominant satellite genome-wide (`satellite_map.py`) |
| `tidk` | ≥0.2.6 | Telomere-repeat windows (`telomere_crosscheck.py`) |
| `nextflow` | ≥23.10 | Pipeline runtime (DSL2), isolated `nextflow` env |

Tool usage is grounded in the offline `docs-corpus/` snapshots (version-matched `--help`/`man` captures, upstream repo docs, and web docs as fallback) — the agent reads these instead of guessing flags.

## 🔄 Update check

Every skill built with [skill-forge](https://github.com/cheahhl814/skill-forge) ships a self-update check that compares the deployed git SHA against this upstream repo via `git fetch` — no GitHub API call, no extra dependencies.

```bash
pixi run update-check
```

| Verdict | Exit | Meaning |
|:--------|:-----|:--------|
| `UP-TO-DATE` | 0 | Local HEAD matches origin/HEAD |
| `LOCAL-AHEAD` | 0 | Unpushed local commits; no action needed |
| `BEHIND-BY-N` | 1 | Upstream is N commits ahead → re-pull/re-sync from this repo |
| `OFFLINE` | 2 | `git fetch` failed; informational only |
| `NO-ORIGIN` | 2 | No `origin` remote configured; informational only |

## 🧩 Depends on

This skill declares hard dependencies on sibling skills (frontmatter `skills:` field — validated by `skill-lint.py` C12; runtime prerequisites live in the separate `requires:` field):

- **read-qc-trimming** — upstream QC and trimming of the nuclear reads used as TideHunter input
- **chartz** — renders the per-scaffold centromere / telomere evidence figures

No reference database is downloaded: all tools consume only user-supplied files, and the only network activity is pixi package installation on first use.

## 📁 Repository layout

```text
centromere-prediction/
├── SKILL.md                    # Master orchestrator (router — start here)
├── README.md                   # This file
├── pixi.toml                   # Pinned tool environment (pixi install) + isolated nextflow env
├── params.json                 # skill-forge render input (tracked)
├── docs-corpus/                # Offline snapshots of upstream tool docs
│   ├── bedtools/  cd-hit/  nextflow/  samtools/
│   └── stringdecomposer/  tidehunter/  trf/
├── preflight/run-preflight/    # Phase 1 sub-skill
├── run/pipeline-run/           # Phase 2 sub-skill
├── qc/results-qc/              # Phase 3 sub-skill
├── debug/run-debug/            # Failure interpreter (signature library)
├── bin/
│   ├── extract_candidate_region.sh   # Densest window per scaffold (window 0 = legacy merged top-10)
│   ├── rank_monomers.py              # cd-hit clusters → top-N by size + ranking table
│   ├── decomp_stats.py               # Per-scaffold dominant-monomer stats
│   ├── satellite_map.py              # BLASTn a monomer genome-wide → arrays
│   ├── telomere_crosscheck.py        # tidk: internal vs capped vs uncapped arrays
│   ├── evidence_tracks.py            # chartz per-scaffold evidence figures
│   └── skill-update-check.py         # Self-update check (pixi run update-check)
├── test-data/                  # Deterministic fixtures + MANIFEST.tsv (+ satellite fixture, MANIFEST_sat.tsv)
├── docs/                       # Workflow metro-map (SVG + PNG, light and dark)
└── runners/nextflow-runner/    # Thin DSL2 runner (main.nf, modules, conf/, test_smoke.py)
```

## 🔒 Hard guarantees

- **Filesystem evidence chain** — every phase emits the artifact the next phase consumes; the boundary between phases is the filesystem, not agent memory.
- **Docs-grounded tool usage** — tool flags come from `docs-corpus/`, never from the agent's memory.
- **Explicit stop points** — ambiguous decisions are surfaced to you as *Evidence + Recommend + Options*, not auto-picked.
- **Reproducible environments** — every tool is pinned in `pixi.toml` and resolved via pixi.
- **Run on real assemblies** — the pipeline has run end-to-end on two real assemblies (173 and 33 scaffolds). `-profile test -stub-run` (empty-signal fixtures) and `-profile test_sat` (a real run on a satellite fixture) are the regression canaries, backed by deterministic fixtures (`test-data/`, fixed seed, sha256 manifest).
- **Honest verdicts** — a negative result from the default 50–250 bp range is provisional until re-checked at `period_max: 2000`, and a satellite array is called a *candidate* centromere until independent evidence (CENP-A ChIP, Hi-C, cytology) supports it.
- **No stale-run contamination** — `publishDir` uses `overwrite: false`; the preflight stale-results check and the run/debug instructions guard against it.

## 🚀 Nextflow runner

For production / HPC / cohort runs, the skill ships a thin Nextflow DSL2 runner at `runners/nextflow-runner/` (resource labels, pinned tool env, trace/report/timeline enabled).

```bash
# Stub smoke test (builds the DAG without executing)
pixi run -e nextflow nextflow run runners/nextflow-runner/main.nf -profile test -stub-run

# Real run on the satellite fixture (exercises every process)
pixi run -e nextflow nextflow run runners/nextflow-runner/main.nf -profile test_sat

# Structural smoke test of the whole skill
python3 runners/nextflow-runner/test_smoke.py
```

The bash recipes in the sub-skills remain the source of truth; the runner is a thin executor.

## Provenance

Built with the [skill-forge](https://github.com/cheahhl814/skill-forge) meta-skill, following the AiX-BIO skill convention (preflight → run → qc → debug evidence chain). Pattern adopted from:

- **BettaMt-agents** — https://github.com/cheahhl814/BettaMt-agents
- **Betta-WGS-agent** — https://github.com/cheahhl814/Betta-WGS-agent
- **bacterial-genome-analysis** — https://github.com/cheahhl814/bacterial-genome-analysis

Worked examples: *Parosphromenus rubrimontis* (500 Mb, 173 scaffolds; TideHunter ~13 min + TRF ~26 min in parallel, under 1 h on 24 cpus / 96 GB) and *Betta coccina* (478 Mb, 33 scaffolds; the run that motivated v1.1.0).

## License

No license file has been published in this repository yet.
