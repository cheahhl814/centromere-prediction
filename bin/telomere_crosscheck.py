#!/usr/bin/env python3
"""Cross-check satellite arrays against telomere repeats (tidk) to separate centromeric from sub-telomeric.

Run inside the skill env so tidk is on PATH:
  pixi run --manifest-path <skill-root>/pixi.toml python3 bin/telomere_crosscheck.py \\
      --assembly asm.fasta --blocks satellite/satellite_main_blocks.tsv --outdir telomere

Outputs: tidk_telomeric_repeat_windows.tsv (the tidk search table; reused if present, or pass --tidk-tsv),
block_vs_telomere.tsv. Call per scaffold:
  internal          block >= --internal-min bp from the nearest end (a sub-telomeric satellite sits at a telomere)
  near end, capped  block within reach of an end that carries telomere repeat
  near end, uncapped  block within reach of an end with NO telomere repeat: the satellite abuts an assembly
                    break (typical of an acrocentric chromosome with an unassembled short arm)
Telomere presence = >= --tel-min repeats in the terminal --terminal-kb of that end. The default motif is the
vertebrate TTAGGG; check the taxon before trusting an "uncapped" call. Stdlib only.
"""
import argparse
import collections
import os
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--blocks", required=True, help="satellite_main_blocks.tsv from satellite_map.py")
    ap.add_argument("--outdir", default="telomere")
    ap.add_argument("--motif", default="TTAGGG")
    ap.add_argument("--window", type=int, default=10000)
    ap.add_argument("--tidk-tsv", help="reuse an existing tidk search table instead of running tidk")
    ap.add_argument("--terminal-kb", type=int, default=50)
    ap.add_argument("--tel-min", type=int, default=100)
    ap.add_argument("--internal-min", type=int, default=2000000)
    ap.add_argument("--min-block-hits", type=int, default=30, help="skip weak blocks")
    ap.add_argument("--min-scaffold", type=int, default=5000000, help="skip scaffolds shorter than this (bp)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    tsv = a.tidk_tsv or os.path.join(a.outdir, "tidk_telomeric_repeat_windows.tsv")  # tidk names it <-o>_telomeric_repeat_windows.tsv
    if not os.path.exists(tsv):
        subprocess.run(["tidk", "search", "-s", a.motif, "-w", str(a.window), "-o", "tidk", "-d", a.outdir,
                        "--extension", "tsv", a.assembly], check=True, stdout=subprocess.DEVNULL)
    cnt = collections.defaultdict(dict)  # window END coordinate -> repeats
    for line in list(open(tsv))[1:]:
        c, win, f, r = line.split("\t")[:4]
        cnt[c][int(win)] = int(f) + int(r)
    tel = lambda c, lo, hi: sum(v for k, v in cnt[c].items() if lo < k <= hi)
    out = ["scaffold\tlen_Mb\tblock_Mb\tdist_left_kb\tdist_right_kb\ttel_left\ttel_right\tcall"]
    for line in list(open(a.blocks))[1:]:
        c, ln, s, e, span, nh, pid, pos = line.rstrip("\n").split("\t")
        ln, s, e = int(ln), int(s), int(e)
        if int(nh) < a.min_block_hits or ln < a.min_scaffold:
            continue
        dl, dr, tk = s, ln - e, a.terminal_kb * 1000
        tl, tr = tel(c, 0, tk), tel(c, ln - tk, ln + a.window)
        side, near, ts = ("left", dl, tl) if dl < dr else ("right", dr, tr)
        if near >= a.internal_min:
            call = "internal"
        elif ts >= a.tel_min:
            call = f"near end, capped ({side})"
        else:
            call = f"near end, uncapped ({side})"
        out.append(f"{c}\t{ln / 1e6:.1f}\t{s / 1e6:.2f}-{e / 1e6:.2f}\t{dl // 1000}\t{dr // 1000}\t{tl}\t{tr}\t{call}")
    path = os.path.join(a.outdir, "block_vs_telomere.tsv")
    open(path, "w").write("\n".join(out) + "\n")
    both = sum(1 for l in out[1:] if int(l.split("\t")[5]) >= a.tel_min and int(l.split("\t")[6]) >= a.tel_min)
    print(f"{len(out) - 1} scaffolds -> {path}; {both} have telomere repeat at both ends", file=sys.stderr)


if __name__ == "__main__":
    main()
