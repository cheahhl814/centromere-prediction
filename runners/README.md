# runners/ — the thin Nextflow DSL2 runner

`nextflow-runner/` is the executable pipeline (5 processes in
`modules/centromere.nf`, kept verbatim from the proven standalone
`centromere-prediction` pipeline — commit `143a5dc` value-channel fix included).
The bash recipes in the sub-skills (`preflight/`, `run/`, `qc/`, `debug/`) are
the source of truth for the procedure; the runner is the bulk executor.

Layout:

| File | Purpose |
|---|---|
| `main.nf` | entry — validates params, creates value/queue channels |
| `modules/centromere.nf` | the 5 processes |
| `nextflow.config` | params, profiles (local/slurm), process resources, trace/report/timeline |
| `conf/test.config` | test profile — stub-run against `test-data/` fixtures |
| `params.example.json` | worked-example params (replace paths before use) |
| `test_smoke.py` | stub-run smoke test (no tool execution) |

Run:

```bash
NXF="nextflow"   # or: pixi exec -e nextflow nextflow
nextflow run nextflow-runner/main.nf -profile local -params-file <your params>
# CI stub run:
nextflow run nextflow-runner/main.nf -profile test -stub-run
```

Tool env: shared one-manifest pixi env at the skill root
(`pixi run --manifest-path <skill-root>/pixi.toml <tool>`); module manifests
point there via `${projectDir}/../../pixi.toml`. Nextflow itself: `pixi run -e
nextflow nextflow ...`.
