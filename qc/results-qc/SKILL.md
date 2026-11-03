---
name: centromere-prediction-results-qc
description: 'Phase 3 of centromere-prediction: interpret results/centromere/ outputs (monomer ranking, TRF density windows, candidate regions, stringDecomposer decomposition, per-scaffold summary), then map the dominant satellite genome-wide, cross-check it against telomere repeats, plot the evidence, and write centromere_report.md with a verdict (good | partial | fragmented | monomer_filter_too_strict | no_centromeric_signal). Use when the centromere run finished, or for checking and interpreting existing centromere results.'
version: 1.1.0
updated: "2026-11-03"
triggers:
  - "check centromere results"
  - "interpret centromere"
  - "centromere report"
  - "map the satellite"
  - "centromere telomere"
---

# Sub-skill: results-qc (Phase 3 — The Report)

## 0. Inputs / Outputs contract

| Direction | Artifact | Notes |
| --- | --- | --- |
| Inputs | `<RUN_DIR>/results/centromere/**` | from pipeline-run (or any existing results dir) |
| Outputs | `<RUN_DIR>/centromere_report.md` with verdict + warnings | consumed by the user / downstream analysis |
| Outputs | `satellite/`, `telomere/`, `charts/` (steps 6–8) | optional but expected before a `good`/`partial` verdict |

Refuse if `results/centromere/` is missing entirely — the run never happened (route to run sub-skill).

## Expected layout

```
results/centromere/
├── centromere_summary.tsv          # scaffold  start  end  n_repeats   (n_repeats = number of TRF entries, NOT copies)
├── regions/candidate_regions.tsv   # chr  start  end  (densest --region_window per scaffold, 0-based start)
├── tidehunter/
│   ├── consensus_clustered.fasta   # all cd-hit clusters (representatives)
│   ├── monomers_top.fasta          # the --max_monomers most abundant clusters (what stringDecomposer used)
│   └── monomer_ranking.tsv         # rank  monomer_id  cluster_size  length  kept
├── trf/
│   ├── all_repeats.bed                    # unfiltered TRF output (period <= --period_max)
│   ├── centromere_candidate_repeats.bed   # period_min–period_max bp, >= min_copies
│   ├── centromere_density_100kb.bed       # repeat density per 100 kb window (score = col 7)
│   └── genome.sizes
└── stringdecomposer/<scaffold>/
    ├── decomposition.tsv                  # region  monomer(' = reverse)  start  end  identity  ...  (v1.1.0+; v1.0.0 held only the log)
    ├── stringdecomposer.log
    └── <scaffold>_centromere.fa           # extracted candidate region
```

`CENTRO_DIR="$RUN_DIR/results/centromere"` below. Results from v1.0.0 have a log, not a table, in `decomposition.tsv`: rerun, or read `decomposition_dir/final_decomposition.tsv` from the Nextflow work dir.

## Procedure — the interpretations

### 1. Monomers (TideHunter + cd-hit-est)

```bash
sort -k3,3nr "$CENTRO_DIR/tidehunter/monomer_ranking.tsv" | head -10   # biggest clusters first
```

The cluster **count** is not diagnostic: long reads give thousands (3,017–3,186 in *B. coccina*) and most are singletons or noise. What matters is whether a few families stand out — a cluster size far above the median, and long monomers. Weak signal if no cluster stands out; mtDNA/organellar repeats pollute untrimmed reads.

### 2. TRF density and candidate regions

```bash
sort -k7,7nr "$CENTRO_DIR/trf/centromere_density_100kb.bed" | head -10
```

Column 7 = repeat bp per window bp (0–1). **TRF reports a clean array as ONE long repeat**, so a window at density 1.0 with 1–3 repeats is a real long array, not an artifact (verified on a 60 kb satellite fixture → a single 59,851 bp entry). A candidate region is the densest window, which may not be the functional centromere.

### 3. stringDecomposer structure

```bash
python3 -I <skill-root>/bin/decomp_stats.py "$CENTRO_DIR" decomp   # -> decomp.scaffold_stats.tsv, decomp.monomer_dominance.tsv
```

- **Dominant monomer:** `top_share` (fraction of decomposed bases on the top monomer) and `top_mean_ident` (identity of calls to that monomer). One monomer ≥ 80% over the region is a high-quality satellite; 25–70% with `top_mean_ident` ≥ 75 is a real but diverged or higher-order array (*B. coccina*: 63% and 77.8% for a 685 bp unit = ~10 × 69 bp); ≤ 15% across 60–100 distinct monomers is background.
- **Do not use `median_ident_all`**: it sits at ~50–60% even where a real satellite exists, because short low-scoring filler calls dominate the count.
- `monomer_dominance.tsv` ranks monomers by how many scaffolds they hold ≥ 25% of; raw "top on N scaffolds" is not used (a weak monomer tops many scaffolds at ~8%).

### 4. Summary

```bash
column -t "$CENTRO_DIR/centromere_summary.tsv" | head -40
```

Columns: `scaffold  start  end  n_repeats` (union span of candidate repeats; `n_repeats` counts TRF entries, so a single clean array can be 1). Healthy: rows for the chromosome-scale scaffolds. Red flags: > 50 rows → fragmented assembly; zero rows → nothing in the period range (go to §Period range).

### 5. Second pass before any negative verdict

If the default 50–250 bp run finds no dominant monomer, **rerun with `--period_max 2000`** (see `run/pipeline-run`) before concluding. TRF stops at `period_max` itself, so the narrow run cannot see longer satellites: in *B. coccina* it missed a 685 bp unit and a ~1.2 kb family that the 2000 bp run found. A negative verdict from the default range is provisional.

### 6. Map the dominant satellite genome-wide

```bash
S=<skill-root>; P="pixi run --manifest-path $S/pixi.toml"
head -3 decomp.monomer_dominance.tsv        # pick the first monomer id (prefix is enough)
$P python3 -I $S/bin/satellite_map.py --assembly <asm.fasta> \
    --monomers "$CENTRO_DIR/tidehunter/consensus_clustered.fasta" --monomer-id <prefix> --outdir satellite
```

Outputs `satellite_blocks.tsv` and `satellite_main_blocks.tsv` (one block per scaffold: span, hit count, mean identity, relative position). A single array of 100–700 kb on most chromosome-scale scaffolds is the centromeric-satellite pattern; scattered short hits are not. Higher identity (85–93%) marks the best-conserved arrays; ~80% is more diverged.

### 7. Telomere cross-check

```bash
$P python3 -I $S/bin/telomere_crosscheck.py --assembly <asm.fasta> --blocks satellite/satellite_main_blocks.tsv --outdir telomere
```

Calls per scaffold: `internal` (≥ 2 Mb from either end: not sub-telomeric); `near end, capped` (block near an end that has telomere repeat: possibly sub-telomeric — inspect); `near end, uncapped` (satellite abuts an assembly break with no telomere there: consistent with a centromere on an acrocentric chromosome whose short arm is unassembled). Default motif TTAGGG (vertebrate): confirm the taxon first.

### 8. Plot the evidence (chartz GenomeTracks)

```bash
python3 -I $S/bin/evidence_tracks.py --density "$CENTRO_DIR/trf/centromere_density_100kb.bed" \
  --hits satellite/satellite_hits.tsv --blocks satellite/satellite_blocks.tsv --calls telomere/block_vs_telomere.tsv \
  --tidk telomere/tidk_telomeric_repeat_windows.tsv --scaffold-lengths satellite/satellite_main_blocks.tsv --outdir charts
```

Open `charts/centromere_telomere_index.html`. Four tracks on shared y-scales: TRF density, satellite coverage, TTAGGG repeats, arrays.

### 9. Report + verdict

Write `centromere_report.md` with: inputs and any deviation (read subsample, capped monomers/window, period range used); monomer ranking highlights; per-scaffold decomposition stats; satellite map and telomere calls; the verdict.

| Verdict | Evidence |
|---|---|
| `good` | one satellite family forms a single array (≥ 100 kb) on ≳ half of the chromosome-scale scaffolds, mostly `internal` or `near end, uncapped` |
| `partial` | satellite evidence on only a few scaffolds (e.g. 2 of ~20): report them, and say what would find the rest (wider window, longer period, other monomers) |
| `fragmented` | > 50 summary rows — candidates spread over contigs; scaffold first |
| `monomer_filter_too_strict` | near-zero repeats at the **wide** period range though biology suggests centromeres |
| `no_centromeric_signal` | no dominant monomer and no satellite array even at `--period_max 2000` — verify taxonomy assumptions |

## What NOT to do

- Do **not** trust `candidate_regions.tsv` as "the centromere" — it is the densest window, not a proven centromere.
- Do **not** call a satellite array a centromere without independent evidence (CENP-A ChIP, Hi-C, cytology): say "candidate".
- Do **not** assume centromeres are mid-chromosome; many species are acrocentric.
- Do **not** conclude "no signal" from the default 50–250 bp run alone (§5).

## Period range

`--period_min/--period_max` (default 50–250; TRF caps at 2000) apply to **both** TideHunter and TRF. A post-filter such as `awk '$7>=pmin && $7<=pmax'` on `centromere_candidate_repeats.bed` can only narrow the range, because TRF itself stops at `period_max`; widening needs a rerun. Guidance: vertebrate/fish centromere satellites are usually 100–2000 bp (use 2000); plants ~150–180 bp or 500+; Drosophila (~5–10 bp) and holocentric organisms (*C. elegans*) are out of scope.

## Signature library

| Signature in QC | Likely cause | Fix |
| --- | --- | --- |
| `decomposition.tsv` has only a few log lines | v1.0.0 published stringDecomposer's log, not the table | rerun with v1.1.0, or use `decomposition_dir/final_decomposition.tsv` from the work dir |
| `median_ident_all` ≈ 55–60% everywhere | filler calls dominate the median, not evidence of "no signal" | use `top_share` + `top_mean_ident` (`decomp_stats.py`) |
| density 1.0 windows with 1–3 repeats | TRF reports a clean array as one repeat | real long arrays; check them, don't discard |
| `monomer_dominance.tsv` top row has ~8% share | weak monomer that tops many scaffolds | read the `scaffolds_top_ge_25pct` column, not `scaffolds_top` |
| `centromere_density_100kb.bed` empty | TRF produced no output (module swallows errors) | rerun TRF standalone, inspect stderr |
| `consensus_clustered.fasta` empty | TideHunter found no tandem repeats | long reads (> 10 kb ideal); check `--period_min/--period_max` |
| `stringdecomposer/<scaffold>/decomposition.tsv` empty + log `skipped: no monomers` | zero monomers (aligner segfaults on an empty set) | fix upstream: reads/period range |
| `stringdecomposer/` missing for a scaffold | no non-zero density | expected, not an error |
| `centromere_summary.tsv` empty | nothing in the period range | widen `--period_max` → `monomer_filter_too_strict` only if still empty |
| tidk table empty / all-zero tel counts | wrong motif for the taxon | check the telomere motif for the clade, rerun with `--motif` |
