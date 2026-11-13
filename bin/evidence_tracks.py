#!/usr/bin/env python3
"""Per-scaffold centromere-vs-telomere evidence figures (chartz GenomeTracks engine) + an index page.

Usage:
  python3 bin/evidence_tracks.py --density results/centromere/trf/centromere_density_100kb.bed \\
      --hits satellite/satellite_hits.tsv --blocks satellite/satellite_blocks.tsv \\
      --calls telomere/block_vs_telomere.tsv --tidk telomere/tidk_telomeric_repeat_windows.tsv \\
      --scaffold-lengths satellite/satellite_main_blocks.tsv --outdir charts

One figure per scaffold listed in --calls, four tracks on shared y-scales so scaffolds compare:
TRF density (0-1), satellite coverage (fraction of each 100 kb window covered by merged BLAST hits,
0-1), TTAGGG repeats per 100 kb (capped at --tel-max), and the merged satellite arrays. Fills the
chartz genetracks template (pinned, zero-dependency SVG; Download SVG/PNG built in). Stdlib only.
--template defaults to ../chartz/templates/genetracks.html next to this skill.
"""
import argparse
import collections
import json
import os
import sys

WIN = 100000


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--density", required=True)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--blocks", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--tidk", required=True)
    ap.add_argument("--scaffold-lengths", required=True, help="a satellite_*blocks.tsv (has scaffold, scaffold_len)")
    ap.add_argument("--outdir", default="charts")
    ap.add_argument("--tel-max", type=int, default=1000)
    ap.add_argument("--min-block-hits", type=int, default=5)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--template", default=os.path.normpath(os.path.join(here, "..", "..", "chartz", "templates", "genetracks.html")))
    a = ap.parse_args()
    tpl = open(a.template).read()
    if tpl.count("/*__SPEC__*/") != 1:
        sys.exit("template must contain exactly one /*__SPEC__*/ slot: " + a.template)
    os.makedirs(a.outdir, exist_ok=True)

    L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in list(open(a.scaffold_lengths))[1:]}
    dens = collections.defaultdict(dict)
    for l in open(a.density):
        f = l.split("\t")
        dens[f[0]][int(f[1])] = float(f[6])
    hits = collections.defaultdict(list)
    for l in list(open(a.hits))[1:]:
        c, s, e, _ = l.split("\t")
        hits[c].append((int(s), int(e)))
    sat = collections.defaultdict(lambda: collections.defaultdict(int))  # bp covered per window, hits merged
    for c, hs in hits.items():
        hs.sort()
        merged = []
        for s, e in hs:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        for s, e in merged:
            for w in range(s // WIN * WIN, e, WIN):
                sat[c][w] += min(e, w + WIN) - max(s, w)
    tel = collections.defaultdict(lambda: collections.defaultdict(int))
    for l in list(open(a.tidk))[1:]:
        c, win, f, r = l.split("\t")[:4]
        tel[c][(int(win) - 1) // WIN * WIN] += int(f) + int(r)  # tidk windows are labelled by their END coordinate
    blocks = collections.defaultdict(list)
    for l in list(open(a.blocks))[1:]:
        c, _, s, e, _, n, _, _ = l.split("\t")
        if int(n) >= a.min_block_hits:
            blocks[c].append((int(s), int(e), int(n)))
    calls = [l.rstrip("\n").split("\t") for l in list(open(a.calls))[1:]]
    cap = lambda v: round(min(v, 1.0), 3)
    names = []
    for row in calls:
        c, block_mb, call = row[0], row[2], row[-1]
        n = L[c]
        starts = range(0, n, WIN)
        main = max(blocks[c], key=lambda b: b[2])
        spec = {
            "title": f"{c} ({n / 1e6:.1f} Mb): centromere vs telomere evidence (100 kb windows; block {block_mb} Mb; {call})",
            "region": {"chrom": c, "start": 0, "end": n},
            "tracks": [
                {"type": "coverage", "label": "TRF density", "height": 60, "color": "#0072B2", "max": 1,
                 "bins": [[s, cap(dens[c].get(s, 0))] for s in starts]},
                {"type": "coverage", "label": "Satellite (BLAST)", "height": 60, "color": "#D55E00", "max": 1,
                 "bins": [[s, cap(sat[c].get(s, 0) / WIN)] for s in starts]},
                {"type": "coverage", "label": "TTAGGG repeats / 100 kb", "height": 60, "color": "#009E73", "max": a.tel_max,
                 "bins": [[s, tel[c].get(s, 0)] for s in starts]},
                {"type": "peaks", "label": "Satellite arrays", "height": 24, "color": "#CC79A7",
                 "items": [dict(start=b[0], end=b[1], **({"label": "main"} if b == main else {})) for b in blocks[c]]}],
            "export": {"width_mm": 183, "dpi": 300, "font_family": "Arial"}}
        text = json.dumps(spec, separators=(",", ":"))
        assert "NaN" not in text and "undefined" not in text
        out = tpl.replace("/*__SPEC__*/", text)
        assert "__SPEC__" not in out
        open(os.path.join(a.outdir, f"centromere_telomere_{c}.html"), "w").write(out)
        names.append(c)
    cards = "\n".join(f'<section><a href="centromere_telomere_{c}.html">{c}</a><iframe src="centromere_telomere_{c}.html" loading="lazy"></iframe></section>' for c in names)
    open(os.path.join(a.outdir, "centromere_telomere_index.html"), "w").write(f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Centromere vs telomere evidence</title>
<style>body{{font-family:system-ui,sans-serif;margin:0;background:#f8f9fa;color:#333}}main{{max-width:1100px;margin:1.5rem auto;padding:0 1rem}}
h1{{font-size:1.1rem}}p{{font-size:.85rem;color:#555;line-height:1.5}}section{{margin:1rem 0}}a{{font-weight:600;font-size:.9rem}}
iframe{{display:block;width:100%;height:440px;border:1px solid #e3e6ea;border-radius:8px;background:#fff;margin-top:.3rem}}</style></head><body><main>
<h1>Centromere vs telomere evidence, {len(names)} scaffolds with a satellite array</h1>
<p>Per 100 kb window: TRF repeat density, coverage by the mapped satellite (BLASTn hits), and TTAGGG telomere repeats (tidk, capped at {a.tel_max}). Pink bars are merged satellite arrays; the main block is labelled. Y-scales are shared across scaffolds. Satellite arrays are candidate centromeres, not proven.</p>
{cards}
</main></body></html>""")
    print(f"{len(names)} figures + index -> {a.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
