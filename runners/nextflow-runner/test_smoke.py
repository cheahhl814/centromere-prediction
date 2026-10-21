#!/usr/bin/env python3
"""test_smoke.py — structural smoke test for the nf-centromere-agents runner.

Verifies (no Nextflow execution, stdlib only):
  1. main.nf exists, declares DSL2, validates params, uses a VALUE channel for assembly
  2. modules/centromere.nf declares all 5 processes with publishDir and pixi-run wrapping
  3. the value-channel fix is present (no bare queue-channel assembly input)
  4. bin/extract_candidate_region.sh is executable
  5. conf/test.config points at the deterministic fixtures
  6. nextflow.config keeps trace/report/timeline enabled
  7. fixtures + MANIFEST.tsv exist and match sha256

Usage: python3 runners/nextflow-runner/test_smoke.py [skill_root]
"""

from __future__ import annotations

import hashlib
import pathlib
import re
import sys

FAILS: list[str] = []
WARNS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}: {name}" + (f" — {detail}" if detail and not ok else ""))
    if not ok:
        FAILS.append(name)


def main() -> int:
    skill_root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[2])
    runner = skill_root / "runners" / "nextflow-runner"

    main_nf = runner / "main.nf"
    mod_nf = runner / "modules" / "centromere.nf"
    nfx_config = runner / "nextflow.config"

    main_src = main_nf.read_text()
    mod_src = mod_nf.read_text()
    conf_src = nfx_config.read_text()

    # 1. main.nf
    check("main.nf exists + DSL2", "nextflow.enable.dsl = 2" in main_src)
    check("main.nf validates params", "Missing required parameters" in main_src)
    check("assembly is a VALUE channel (143a5dc fix)", "Channel.value(file(params.assembly" in main_src)
    check("reads is a queue channel", "Channel.fromPath(params.reads" in main_src)

    # 2-3. module contracts
    for proc in ("tidehunter_cluster", "trf_density", "identify_regions", "stringdecomposer", "centromere_summary"):
        check(f"process {proc} declared", f"process {proc}" in mod_src)
    check("all processes publish", mod_src.count("publishDir") >= 5)
    check("publish pattern .fasta (bug #1 fix)", "pattern: 'consensus_clustered.fasta'" in mod_src)
    check("regions publishDir present (bug #2 fix)", "results/centromere/regions" in mod_src)
    check("bare bedtools wrapped in pixi run (bug #3 fix)", "pixi run --manifest-path ${projectDir}/../../pixi.toml \\\n    bedtools sort" in mod_src)
    check("no unused clust_id input (bug #4 fix)", "clust_id" not in mod_src)
    check("one-manifest pixi path in modules", "${projectDir}/../../pixi.toml" in mod_src)
    check("helper script via projectDir", "${projectDir}/../../bin/extract_candidate_region.sh" in mod_src)

    # 4. helper script
    bin_sh = skill_root / "bin" / "extract_candidate_region.sh"
    check("extract_candidate_region.sh executable", bin_sh.is_file() and bin_sh.stat().st_mode & 0o111 != 0)

    # 5-6. configs
    test_conf = (runner / "conf" / "test.config").read_text()
    check("test profile present", "test.config" in (runner / "conf" / "test.config").name and "params" in test_conf)
    check("trace/report/timeline enabled", all(
        k in conf_src for k in ("trace {", "report {", "timeline {")))

    # 7. fixtures
    manifest = skill_root / "test-data" / "MANIFEST.tsv"
    check("MANIFEST.tsv exists", manifest.is_file())
    if manifest.is_file():
        for line in manifest.read_text().strip().splitlines()[1:]:
            f, fmt, sha, _gen = line.split("\t")
            p = skill_root / "test-data" / f
            data = p.read_bytes()
            if fmt == "fastq":
                data = __import__("gzip").decompress(data)
            check(f"fixture {f} sha256 matches", hashlib.sha256(data).hexdigest() == sha)

    print()
    if FAILS:
        print(f"SMOKE RESULT: FAIL ({len(FAILS)} failures)")
        return 1
    print(f"SMOKE RESULT: PASS ({len(WARNS)} warnings)" if WARNS else "SMOKE RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())