#!/usr/bin/env python3
"""Map one candidate satellite monomer across the assembly (BLASTn) and merge hits into arrays.

Run inside the skill env so blastn/makeblastdb are on PATH:
  pixi run --manifest-path <skill-root>/pixi.toml python3 bin/satellite_map.py \\
      --assembly asm.fasta --monomers results/centromere/tidehunter/consensus_clustered.fasta \\
      --monomer-id 0829ee6a --outdir satellite

--monomer-id is an id PREFIX of a record in --monomers (take it from decomp_stats.py's
monomer_dominance.tsv). Outputs in --outdir: satellite_hits.tsv (filtered BLAST hits),
satellite_blocks.tsv (hits merged within --merge-gap, >= --min-block-hits hits) and
satellite_main_blocks.tsv (the block with most hits per scaffold). Stdlib only.
"""
import argparse
import collections
import os
import subprocess
import sys


def fasta_lengths(path):
    fai = path + ".fai"
    if os.path.exists(fai):
        return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = 0
        elif name:
            out[name] += len(line.strip())
    return out


def pick_monomer(path, prefix):
    name, buf, hits = None, [], []
    for line in open(path):
        if line.startswith(">"):
            if name and name.startswith(prefix):
                hits.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name and name.startswith(prefix):
        hits.append((name, "".join(buf)))
    if len(hits) != 1:
        sys.exit(f"--monomer-id {prefix!r} matched {len(hits)} records (need exactly 1)")
    return hits[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--monomers", required=True)
    ap.add_argument("--monomer-id", required=True)
    ap.add_argument("--outdir", default="satellite")
    ap.add_argument("--evalue", default="1e-10")
    ap.add_argument("--pid", default="70", help="minimum percent identity")
    ap.add_argument("--min-hit", type=int, default=300, help="drop hits shorter than this (bp); removes sub-unit matches")
    ap.add_argument("--merge-gap", type=int, default=20000)
    ap.add_argument("--min-block-hits", type=int, default=5)
    ap.add_argument("--threads", default="16")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    name, seq = pick_monomer(a.monomers, a.monomer_id)
    q = os.path.join(a.outdir, "monomer.fa")
    open(q, "w").write(f">{name}\n{seq}\n")
    db = os.path.join(a.outdir, "genomedb")
    if not os.path.exists(db + ".nsq") and not os.path.exists(db + ".00.nsq"):
        subprocess.run(["makeblastdb", "-in", a.assembly, "-dbtype", "nucl", "-out", db], check=True, stdout=subprocess.DEVNULL)
    raw = subprocess.run(
        ["blastn", "-query", q, "-db", db, "-task", "blastn", "-evalue", a.evalue, "-perc_identity", a.pid,
         "-num_threads", a.threads, "-max_target_seqs", "100000",
         "-outfmt", "6 sseqid pident length sstart send"], check=True, capture_output=True, text=True).stdout
    hits = collections.defaultdict(list)
    for line in raw.splitlines():
        sid, pid, ln, s, e = line.split("\t")
        if int(ln) >= a.min_hit:
            hits[sid].append((min(int(s), int(e)) - 1, max(int(s), int(e)), float(pid)))
    L = fasta_lengths(a.assembly)
    with open(os.path.join(a.outdir, "satellite_hits.tsv"), "w") as fh:
        fh.write("scaffold\tstart\tend\tpident\n")
        for sc in sorted(hits):
            for s, e, p in sorted(hits[sc]):
                fh.write(f"{sc}\t{s}\t{e}\t{p}\n")
    blocks = []
    for sc, hs in hits.items():
        hs.sort()
        cur = None
        for s, e, p in hs:
            if cur and s - cur[1] <= a.merge_gap:
                cur[1] = max(cur[1], e); cur[2] += 1; cur[3] += p
            else:
                if cur: blocks.append(cur)
                cur = [s, e, 1, p, sc]
        if cur: blocks.append(cur)
    blocks = [[b[4], b[0], b[1], b[2], b[3] / b[2]] for b in blocks if b[2] >= a.min_block_hits]
    hdr = "scaffold\tscaffold_len\tstart\tend\tspan_bp\tn_hits\tmean_pident\tpos_frac\n"
    row = lambda b: f"{b[0]}\t{L[b[0]]}\t{b[1]}\t{b[2]}\t{b[2]-b[1]}\t{b[3]}\t{b[4]:.1f}\t{(b[1]+b[2])/2/L[b[0]]:.2f}\n"
    blocks.sort(key=lambda b: (b[0], b[1]))
    with open(os.path.join(a.outdir, "satellite_blocks.tsv"), "w") as fh:
        fh.write(hdr); fh.writelines(row(b) for b in blocks)
    best = {}
    for b in blocks:
        if b[0] not in best or b[3] > best[b[0]][3]:
            best[b[0]] = b
    with open(os.path.join(a.outdir, "satellite_main_blocks.tsv"), "w") as fh:
        fh.write(hdr); fh.writelines(row(b) for b in sorted(best.values(), key=lambda b: -b[3]))
    print(f"{name}: {sum(len(v) for v in hits.values())} hits >= {a.min_hit} bp on {len(hits)} scaffolds; "
          f"{len(blocks)} blocks; {len(best)} scaffolds with a main block", file=sys.stderr)


if __name__ == "__main__":
    main()
