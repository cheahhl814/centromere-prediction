#!/usr/bin/env python3
"""Rank cd-hit-est clusters by member count and keep the top N representatives.

stringDecomposer aligns a region against every monomer it is given, so memory and run time
grow with (region length x monomer count). Ranking by cluster size keeps the most abundant
tandem-repeat families and drops the long tail of singletons.

Usage: rank_monomers.py <clustered.fasta> <clustered.fasta.clstr> <top_n> <out_top.fasta> <out_ranking.tsv>

Writes the top-N representatives to <out_top.fasta> (largest cluster first) and a ranking of
ALL clusters to <out_ranking.tsv> (rank, id, cluster_size, length, kept). top_n <= 0 keeps all.
Stdlib only.
"""
import re
import sys


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def read_clusters(path):
    """cd-hit .clstr -> list of (member_count, representative_id)."""
    clusters, count, rep = [], 0, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                if rep is not None:
                    clusters.append((count, rep))
                count, rep = 0, None
            elif line.strip():
                count += 1
                if line.rstrip().endswith("*"):
                    rep = re.search(r">(\S+?)\.\.\.", line).group(1)
        if rep is not None:
            clusters.append((count, rep))
    return clusters


def main(argv):
    if len(argv) != 6:
        sys.exit(__doc__)
    fasta, clstr, top_n, out_fa, out_tsv = argv[1], argv[2], int(argv[3]), argv[4], argv[5]
    seqs = read_fasta(fasta)
    ranked = sorted(read_clusters(clstr), key=lambda c: (-c[0], c[1]))
    keep = ranked if top_n <= 0 else ranked[:top_n]
    kept_ids = {rep for _, rep in keep}
    with open(out_fa, "w") as fa:
        for _, rep in keep:
            fa.write(f">{rep}\n{seqs[rep]}\n")
    with open(out_tsv, "w") as tsv:
        tsv.write("rank\tmonomer_id\tcluster_size\tlength\tkept\n")
        for i, (n, rep) in enumerate(ranked, 1):
            tsv.write(f"{i}\t{rep}\t{n}\t{len(seqs[rep])}\t{int(rep in kept_ids)}\n")
    print(f"clusters={len(ranked)} kept={len(keep)}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)
