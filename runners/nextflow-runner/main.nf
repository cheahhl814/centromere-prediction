#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Workflow parameters
params.reads = ""
// Sequencing reads (fastq/fasta) — nuclear or trimmed reads preferred
params.assembly = ""
// Assembly FASTA (scaffolded or polished)

// Validate required parameters
if (!params.reads || !params.assembly) {
    error "Missing required parameters: provide --reads and --assembly"
}

// Modules
include { CENTROMERE } from "./modules/centromere.nf"

// Workflow
workflow {
    reads_ch    = Channel.fromPath(params.reads, checkIfExists: true)
    // Assembly as a VALUE channel: it is consumed by multiple processes and
    // must be reusable for every stringdecomposer task (queue channels with a
    // single value would limit the process to a single execution).
    assembly_ch = Channel.value(file(params.assembly, checkIfExists: true))

    CENTROMERE(reads_ch, assembly_ch)
}
