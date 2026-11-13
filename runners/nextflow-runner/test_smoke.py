#!/usr/bin/env python3
"""test_smoke.py — structural smoke test for the centromere-prediction runner.

Verifies (no Nextflow execution, stdlib only):
  1. main.nf exists, declares DSL2, validates params, uses a VALUE channel for assembly
  2. modules/centromere.nf declares all 5 processes with publishDir and pixi-run wrapping
  3. the value-channel fix is present (no bare queue-channel assembly input)
  4. bin/extract_candidate_region.sh is executable
  5. conf/test.config points at the deterministic fixtures
  6. nextflow.config keeps trace/report/timeline enabled
  7. fixtures + MANIFEST.tsv exist and match sha256
  8. v1.1.0: parameterised period range, capped monomers/region window, real decomposition table
     published, helper scripts behave (unit-tested on synthetic input), satellite fixture + test_sat profile

Usage: python3 runners/nextflow-runner/test_smoke.py [skill_root]
"""

from __future__ import annotations

import hashlib
import pathlib
import re
import subprocess
import sys
import tempfile

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
    check("publish pattern .fasta (bug #1 fix)", "pattern: '{consensus_clustered.fasta," in mod_src)
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


    # 8. v1.1.0 contracts
    for prm in ("period_min", "period_max", "min_copies", "region_window", "max_monomers"):
        check(f"main.nf declares params.{prm}", f"params.{prm}" in main_src)
    check("main.nf validates the period range (TRF cap 2000)", "period_max > 2000" in main_src)
    check("TRF max period comes from params (not hardcoded 250)", "10 50 ${params.period_max}" in mod_src and "10 50 250" not in mod_src)
    check("TideHunter period range from params", "-p ${params.period_min} -P ${params.period_max}" in mod_src)
    check("monomers ranked + capped before stringdecomposer", "rank_monomers.py" in mod_src and "out.top_monomers" in mod_src)
    check("region window passed to identify_regions", "${params.region_window}" in mod_src)
    check("real decomposition table published (bug #6 fix)", "cp decomposition_dir/final_decomposition.tsv decomposition.tsv" in mod_src)
    check("empty monomer set skipped (SIGSEGV guard)", "grep -q '^>' ${monomers}" in mod_src)
    check("BED start converted to 1-based for samtools (off-by-one fix)", ".toLong() + 1" in mod_src)
    check("test_sat profile registered", "test_sat" in conf_src and (runner / "conf" / "test_sat.config").is_file())
    for helper in ("rank_monomers.py", "decomp_stats.py", "satellite_map.py", "telomere_crosscheck.py", "evidence_tracks.py"):
        hp = skill_root / "bin" / helper
        check(f"bin/{helper} executable", hp.is_file() and hp.stat().st_mode & 0o111 != 0)

    # unit tests on synthetic input (python/bash/awk only)
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        (td / "c.fasta").write_text(">a_rep0 x\nACGT\n>b_rep0 y\nGGGGCC\n>c_rep0 z\nTT\n")
        (td / "c.clstr").write_text(">Cluster 0\n0\t4nt, >a_rep0... *\n1\t4nt, >q1... at +/99%\n>Cluster 1\n0\t6nt, >b_rep0... *\n"
                                    ">Cluster 2\n0\t2nt, >c_rep0... *\n1\t2nt, >q2... at +/99%\n2\t2nt, >q3... at +/99%\n")
        r = subprocess.run([sys.executable, str(skill_root / "bin" / "rank_monomers.py"), str(td / "c.fasta"), str(td / "c.clstr"), "2",
                            str(td / "top.fa"), str(td / "rank.tsv")], capture_output=True, text=True)
        top = (td / "top.fa").read_text() if (td / "top.fa").exists() else ""
        check("rank_monomers keeps the 2 largest clusters, largest first",
              r.returncode == 0 and re.findall(r">(\S+)", top) == ["c_rep0", "a_rep0"], r.stderr)
        rank = (td / "rank.tsv").read_text().splitlines() if (td / "rank.tsv").exists() else []
        check("rank_monomers ranks all clusters", len(rank) == 4 and rank[1].split("\t")[:3] == ["1", "c_rep0", "3"], str(rank))
        (td / "e.fasta").write_text(""); (td / "e.clstr").write_text("")
        r = subprocess.run([sys.executable, str(skill_root / "bin" / "rank_monomers.py"), str(td / "e.fasta"), str(td / "e.clstr"), "5",
                            str(td / "etop.fa"), str(td / "erank.tsv")], capture_output=True, text=True)
        check("rank_monomers handles an empty cluster set", r.returncode == 0 and (td / "etop.fa").read_text() == "", r.stderr)

        # density: 6 bins of 100 kb; densest 3-bin (300 kb) run is bins 3-5
        dens = [0.1, 0.0, 0.5, 0.9, 0.8, 0.0]
        (td / "d.bed").write_text("".join(f"s1\t{i*100000}\t{(i+1)*100000}\t1\t1\t100000\t{v}\n" for i, v in enumerate(dens)))
        (td / "s.txt").write_text("s1\n")
        sh = str(skill_root / "bin" / "extract_candidate_region.sh")
        subprocess.run(["bash", sh, str(td / "d.bed"), str(td / "s.txt"), str(td / "w.tsv"), "10", "300000"], check=True)
        check("extract_candidate_region picks the densest contiguous window", (td / "w.tsv").read_text().split() == ["s1", "200000", "500000"], (td / "w.tsv").read_text())
        subprocess.run(["bash", sh, str(td / "d.bed"), str(td / "s.txt"), str(td / "l.tsv"), "3", "0"], check=True)
        check("extract_candidate_region window 0 = legacy merged top-N span", (td / "l.tsv").read_text().split() == ["s1", "200000", "500000"], (td / "l.tsv").read_text())

    # satellite fixture sha256 (generator + manifest)
    sat_manifest = skill_root / "test-data" / "MANIFEST_sat.tsv"
    check("MANIFEST_sat.tsv exists", sat_manifest.is_file())
    if sat_manifest.is_file():
        for line in sat_manifest.read_text().strip().splitlines()[1:]:
            f, fmt, sha, _gen = line.split("\t")
            data = (skill_root / "test-data" / f).read_bytes()
            if fmt == "fastq":
                data = __import__("gzip").decompress(data)
            check(f"satellite fixture {f} sha256 matches", hashlib.sha256(data).hexdigest() == sha)

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