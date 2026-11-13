#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Workflow parameters
params.reads = ""
// Sequencing reads (fastq/fasta) — nuclear or trimmed reads preferred
params.assembly = ""
// Assembly FASTA (scaffolded or polished)
params.period_min    = 50       // monomer period range, bp (TRF itself stops at period_max, max 2000)
params.period_max    = 250
params.min_copies    = 5        // TRF candidate filter: minimum copy number
params.region_window = 1000000  // bp; decompose only the densest window of this size per scaffold (0 = legacy merged span)
params.max_monomers  = 100      // keep the N most abundant cd-hit clusters for stringDecomposer (0 = all)

// Modules
include { CENTROMERE } from "./modules/centromere.nf"

// Workflow
workflow {
    // Nextflow ≥26 forbids bare statements at script top level — the params
    // validation lives inside the workflow (strict-DSL compliant).
    if (!params.reads || !params.assembly) {
        error "Missing required parameters: provide --reads and --assembly"
    }
    if (params.period_min < 1 || params.period_max <= params.period_min || params.period_max > 2000) {
        error "Invalid period range ${params.period_min}-${params.period_max}: need 1 <= period_min < period_max <= 2000 (TRF limit)"
    }
    reads_ch    = Channel.fromPath(params.reads, checkIfExists: true)
    // Assembly as a VALUE channel: it is consumed by multiple processes and
    // must be reusable for every stringdecomposer task (queue channels with a
    // single value would limit the process to a single execution).
    assembly_ch = Channel.value(file(params.assembly, checkIfExists: true))

    CENTROMERE(reads_ch, assembly_ch)
}
