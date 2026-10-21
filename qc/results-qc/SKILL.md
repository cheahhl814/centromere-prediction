---
name: centromere-prediction-results-qc
description: 'Phase 3 of centromere-prediction: interpret results/centromere/ outputs (monomer counts and lengths, density windows, candidate regions, stringDecomposer decomposition structure, per-scaffold summary) with healthy-range targets and red-flag diagnostics, then write centromere_report.md with a verdict (good | fragmented | monomer_filter_too_strict | no_centromeric_signal). Use when the centromere run finished, or for checking and interpreting existing centromere results.'
version: 1.0.0
updated: "2026-10-21"
triggers:
  - "check centromere results"
  - "interpret centromere"
  - "centromere report"
---

# Sub-skill: results-qc (Phase 3 — The Report)

## 0. Inputs / Outputs contract

| Direction | Artifact | Notes |
| --- | --- | --- |
| Inputs | `<RUN_DIR>/results/centromere/**` | from pipeline-run (or any existing results dir) |
| Outputs | `<RUN_DIR>/centromere_report.md` with verdict + warnings | consumed by the user / downstream analysis |

Refuse if `results/centromere/` is missing entirely — the run never happened (route to run sub-skill).

## Expected layout

```
results/centromere/
├── centromere_summary.tsv          # scaffold  start  end  n_repeats
├── regions/candidate_regions.tsv   # chr  start  end  (top-10 density windows per scaffold, merged)
├── tidehunter/consensus_clustered.fasta   # clustered tandem-repeat monomers
├── trf/
│   ├── centromere_candidate_repeats.bed   # TRF repeats, period 50–250 bp, ≥5 copies
│   ├── all_repeats.bed                    # unfiltered TRF output
│   ├── centromere_density_100kb.bed       # repeat density per 100 kb window (score = col 7)
│   └── genome.sizes
└── stringdecomposer/<scaffold>/
    ├── decomposition.tsv                  # position  monomer_id  orientation  score
    └── <scaffold>_centromere.fa           # extracted candidate region
```

## Procedure — the five interpretations

All commands run with `CENTRO_DIR="$RUN_DIR/results/centromere"`.

### 1. Monomers (TideHunter + cd-hit-est)

```bash
MONOMERS="$CENTRO_DIR/tidehunter/consensus_clustered.fasta"
awk '/^>/{if(seq)print name,length(seq); name=$0; seq=""; next} {seq=seq$0} END{if(seq)print name,length(seq)}' \
    "$MONOMERS" | sort -k2 -nr | head -10
```

**Healthy target:** 50–500 monomers of length 50–250 bp.

- Only a few monomers → cd-hit clustering too aggressive (`-c 0.8`) or weak signal
- Thousands of monomers → reads contained many non-centromeric tandem repeats (e.g. mtDNA); prefer nuclear reads

### 2. TRF density and candidate regions

```bash
DENSITY="$CENTRO_DIR/trf/centromere_density_100kb.bed"
awk -F'\t' '{print $0"\t"($7+0)}' "$DENSITY" | sort -k7,7nr | head -10
```

Column 7 = bp-of-repeat per bp-of-window, typically 0.0–0.5. `extract_candidate_region.sh` merges the top-10 windows per scaffold into one span → `candidate_regions.tsv`.

### 3. stringDecomposer structure

```bash
ls "$CENTRO_DIR/stringdecomposer/"
cat "$CENTRO_DIR/stringdecomposer/<scaffold>/decomposition.tsv"
```

Columns: `position  monomer_id  orientation  score`.

- **High-quality centromere:** one dominant monomer covering >80% of the region, periodic spacing matching the monomer length, high HOR structure in vertebrates
- **Weak/absent:** many different monomers, random spacing

### 4. Summary

```bash
column -t "$CENTRO_DIR/centromere_summary.tsv" | head -25
```

Columns: `scaffold  start  end  n_repeats` (union span of candidate repeats; n_repeats = count in period 50–250 bp).

**Healthy target** (typical fish, ~24 chromosomes): ~24 rows, one per chromosome-scale scaffold, each with n_repeats in the hundreds to thousands.

**Red flags:**

- More than 50 rows → assembly fragmented; candidates spread over contigs
- One scaffold with 10× more repeats → possible collapsed centromeres / misassembly
- Zero rows → no 50–250 bp tandem repeats; motif may differ (see period customization below)

### 5. Report + verdict

Write `centromere_report.md` with: n_monomers + length range; n_scaffolds with a candidate vs expected (≈ chromosome number); dominant-monomer evidence from stringdecomposer; and a verdict:

| Verdict | Evidence |
|---|---|
| `good` | ≈ chromosome-number scaffolds with candidates, hundreds–thousands of repeats each, dominant monomer ≥80% coverage in decompositions |
| `fragmented` | >50 summary rows — candidates spread over contigs; assemble/scaffold first |
| `monomer_filter_too_strict` | near-zero monomers/repeats but biology suggests centromeres exist — loosen period range |
| `no_centromeric_signal` | genuinely no 50–250 bp tandem-repeat signal — verify taxonomy assumptions (§period range) before re-running |

## What NOT to do

- Do **not** trust `candidate_regions.tsv` as "the centromere" — it's the highest-density window, which may not be the functional centromere.
- Do **not** call a single dominant monomer a centromere without genomic context (e.g. CENP-A/CenH3 ChIP-seq).
- Do **not** assume centromeres are mid-chromosome — many species have sub-telomeric centromeres.

## Period range customization

The TRF filter (period 50–250 bp, copies ≥5) is hardcoded in `trf_density` and targets vertebrate centromeres.

- **Plants** (e.g. *Arabidopsis*): monomers 150–180 bp (`AthCEN*`) or 500+ bp — widen the upper bound
- **Drosophila**: monomers ~5–10 bp — too small for this filter; needs a separate analysis
- **Holocentric organisms** (e.g. *C. elegans*): no localized centromere — this analysis is uninformative

Post-process without forking the module:

```bash
awk -v pmin=100 -v pmax=500 '$7+0 >= pmin && $7+0 <= pmax' \
    "$CENTRO_DIR/trf/centromere_candidate_repeats.bed" \
    > "$CENTRO_DIR/trf/centromere_candidate_repeats_filtered.bed"
```

## Signature library

| Signature in QC | Likely cause | Fix |
| --- | --- | --- |
| `centromere_density_100kb.bed` empty | TRF produced no output; `\|\| true` in the module swallowed the error | Re-run TRF manually per the run sub-skill; inspect `all_repeats.bed`. |
| `consensus_clustered.fasta` empty | TideHunter found no tandem repeats | Verify reads are long (>10 kb ideal); check `-p 50 -P 250`. |
| `stringdecomposer/` missing for a scaffold | scaffold had no non-zero density | Expected for some scaffolds — not an error. |
| Only one scaffold dir present but many candidates in `candidate_regions.tsv` | pre-`143a5dc` queue-channel bug left stale results | Trash stale outputs; re-run the fixed runner fresh. |
| `centromere_summary.tsv` empty | no repeats in the 50–250 bp period | Loosen the period range, or the organism isn't tandem-repeat-rich → verdict `monomer_filter_too_strict`. |
| density scores all 0.000 but `all_repeats.bed` is populated | repeats filtered out by period/copies in `centromere_candidate_repeats.bed` | Inspect the period distribution: `awk -F'\t' '{print int($7)}' all_repeats.bed \| sort \| uniq -c \| sort -rn \| head`. |