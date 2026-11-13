#!/usr/bin/env python3
"""Per-scaffold stringDecomposer statistics + which monomer dominates across scaffolds.

Usage: decomp_stats.py <results/centromere> [out_prefix]

Reads <results>/stringdecomposer/<scaffold>/decomposition.tsv (v1.1.0+; falls back to
decomposition_dir/final_decomposition.tsv). Columns of that table: region, monomer
(trailing ' = reverse strand), start, end, identity(%), ...

Writes <out_prefix>.scaffold_stats.tsv and <out_prefix>.monomer_dominance.tsv (default prefix
decomp_stats). Stdlib only.

Interpretation notes
- Judge a scaffold by top_share (fraction of decomposed bases assigned to its top monomer) and
  top_mean_ident (mean identity of calls to THAT monomer). median_ident_all is dominated by
  short low-scoring filler calls (~50-60%) even where a real satellite exists, so it is not
  evidence of "no signal".
- monomer_dominance.tsv ranks monomers by the number of scaffolds where they hold >= MIN_SHARE
  (0.25) of the decomposed bases, then by the summed top share. A monomer that is merely "top" on
  many scaffolds with a ~5-10% share is generic background, so raw top-counts are not used.
  Pass the first row's id prefix to satellite_map.py.
"""
import collections
import glob
import os
import re
import statistics as st
import sys


MIN_SHARE = 0.25


def load(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            rows.append((p[0], p[1].rstrip("'"), int(p[2]), int(p[3]), float(p[4])))
    return rows


def main(argv):
    if len(argv) not in (2, 3):
        sys.exit(__doc__)
    res, prefix = argv[1], (argv[2] if len(argv) == 3 else "decomp_stats")
    dominance = collections.defaultdict(lambda: [0, 0, 0.0])  # monomer -> [n_top, n_top_ge_min_share, sum_top_share]
    out = ["scaffold\tregion_bp\tcalls\tcovered_frac\tn_distinct\ttop_monomer\ttop_share\ttop3_share\ttop_mean_ident\tmedian_ident_all\tmedian_call_len"]
    for d in sorted(glob.glob(os.path.join(res, "stringdecomposer", "*"))):
        scaffold = os.path.basename(d)
        f = os.path.join(d, "decomposition.tsv")
        if not (os.path.exists(f) and os.path.getsize(f) > 0):  # v1.0.0 published only the log here
            f = os.path.join(d, "decomposition_dir", "final_decomposition.tsv")
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            continue
        rows = load(f)
        if not rows:
            continue
        m = re.search(r":(\d+)-(\d+)$", rows[0][0])
        region = int(m.group(2)) - int(m.group(1)) + 1 if m else max(r[3] for r in rows) + 1
        bases, idents = collections.Counter(), collections.defaultdict(list)
        for _, mono, s, e, ident in rows:
            bases[mono] += e - s + 1
            idents[mono].append(ident)
        tot = sum(bases.values())
        top = bases.most_common(3)
        share = top[0][1] / tot
        dominance[top[0][0]][0] += 1
        dominance[top[0][0]][1] += share >= MIN_SHARE
        dominance[top[0][0]][2] += share
        out.append("\t".join(map(str, (
            scaffold, region, len(rows), f"{tot / region:.2f}", len(bases), top[0][0][:40], f"{share:.2f}",
            f"{sum(b for _, b in top) / tot:.2f}", f"{st.mean(idents[top[0][0]]):.1f}",
            f"{st.median(r[4] for r in rows):.1f}", st.median(r[3] - r[2] + 1 for r in rows)))))
    with open(prefix + ".scaffold_stats.tsv", "w") as fh:
        fh.write("\n".join(out) + "\n")
    with open(prefix + ".monomer_dominance.tsv", "w") as fh:
        fh.write(f"monomer\tscaffolds_top\tscaffolds_top_ge_{int(MIN_SHARE * 100)}pct\tmean_top_share\n")
        for mono, (n, n_hi, tot_share) in sorted(dominance.items(), key=lambda kv: (-kv[1][1], -kv[1][2])):
            fh.write(f"{mono}\t{n}\t{n_hi}\t{tot_share / n:.2f}\n")
    print(f"{len(out) - 1} scaffolds -> {prefix}.scaffold_stats.tsv, {prefix}.monomer_dominance.tsv", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)
