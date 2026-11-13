#!/usr/bin/env python3
"""centromere-prediction satellite fixture (v1.1.0) — exercises the monomer path.

The v1.0.0 fixtures (make_fixtures.py) contain no tandem repeats, so TideHunter finds
no monomers and the ranking / decomposition steps never run on real data. This
generator writes a tiny assembly with a 171 bp satellite array plus reads drawn from it.

Writes sat_ref.fasta, sat_reads.fastq.gz and MANIFEST_sat.tsv into <outdir>.
Stdlib only, fixed seed, gzip with mtime=0; sha256 is over the DECOMPRESSED bytes.

Usage: python3 make_satellite_fixture.py <outdir>
"""
import gzip
import hashlib
import io
import random
import sys
from pathlib import Path

SEED = 20261103
MONO_LEN = 171


def mutate(rng, seq, rate):
    return "".join(rng.choice("ACGT") if rng.random() < rate else b for b in seq)


def build():
    rng = random.Random(SEED)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    mono = rand(MONO_LEN)
    array = lambda copies: "".join(mutate(rng, mono, 0.06) for _ in range(copies))
    scaffolds = {
        "sat1": (rand(150000), array(350), rand(90000)),   # 60 kb array, density window 100-200 kb
        "sat2": (rand(120000), array(180), rand(50000)),   # 31 kb array
    }
    ref = {name: "".join(parts) for name, parts in scaffolds.items()}
    reads, i = [], 0
    for name, (left, arr, _right) in scaffolds.items():
        for _ in range(45):
            n = rng.randint(4000, 7000)
            start = len(left) + rng.randint(0, max(0, len(arr) - n))
            reads.append(mutate(rng, ref[name][start:start + n], 0.03))
    fa = "".join(f">{n}\n{s}\n" for n, s in ref.items()).encode()
    fq = "".join(f"@r{k}\n{r}\n+\n{'I' * len(r)}\n" for k, r in enumerate(reads, 1)).encode()
    return fa, fq


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    out = Path(argv[1]); out.mkdir(parents=True, exist_ok=True)
    fa, fq = build()
    (out / "sat_ref.fasta").write_bytes(fa)
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, mtime=0) as fh:
        fh.write(fq)
    (out / "sat_reads.fastq.gz").write_bytes(buf.getvalue())
    rows = ["file\tformat\tsha256\tgenerated-by"]
    for f, fmt, data in (("sat_ref.fasta", "fasta", fa), ("sat_reads.fastq.gz", "fastq", fq)):
        rows.append(f"{f}\t{fmt}\t{hashlib.sha256(data).hexdigest()}\tmake_satellite_fixture.py seed={SEED}")
        print(f"OK: {f} ({len(data)} bytes decompressed)")
    (out / "MANIFEST_sat.tsv").write_text("\n".join(rows) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
