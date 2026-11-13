process tidehunter_cluster {
  tag "Discover tandem repeat monomers from reads"
  label 'centromere'
  publishDir "${launchDir}/results/centromere/tidehunter", mode: 'copy', overwrite: false, pattern: '{consensus_clustered.fasta,monomers_top.fasta,monomer_ranking.tsv}'

  input:
  path reads

  output:
  path "consensus_clustered.fasta", emit: monomers
  path "monomers_top.fasta", emit: top_monomers
  path "monomer_ranking.tsv", emit: ranking

  script:
  """
  env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml \\
    TideHunter -p ${params.period_min} -P ${params.period_max} -c 5 -t ${task.cpus} ${reads} > consensus_raw.fasta

  env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml \\
    cd-hit-est -i consensus_raw.fasta -o consensus_clustered.fasta -c 0.8 -n 5 -T ${task.cpus} -M 8000 -d 0 -sc 1 -sf 1

  # stringDecomposer cost grows with (region length x monomer count): keep the most abundant clusters
  env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml \\
    python3 ${projectDir}/../../bin/rank_monomers.py consensus_clustered.fasta consensus_clustered.fasta.clstr ${params.max_monomers} monomers_top.fasta monomer_ranking.tsv
  """
}

process trf_density {
  tag "Map tandem repeats and compute density"
  label 'centromere'
  publishDir "${launchDir}/results/centromere/trf", mode: 'copy', overwrite: false

  input:
  path assembly

  output:
  path "all_repeats.bed", emit: all_repeats_bed
  path "centromere_candidate_repeats.bed", emit: repeats_bed
  path "centromere_density_100kb.bed", emit: density_bed
  path "genome.sizes", emit: sizes

  script:
  """
  env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml \
    trf ${assembly} 2 7 7 80 10 50 ${params.period_max} -f -d -m -ngs -h > assembly.trf.dat 2>/dev/null || true

  awk '/^@/ { seqname = substr(\$1, 2); next }
  NF >= 14 {
      start  = \$1 - 1
      end    = \$2
      period = \$3
      copies = \$4
      score  = \$8
      printf "%s\\t%d\\t%d\\t%s_%d_%s\\t%s\\t.\\t%s\\t%s\\n",
          seqname, start, end, "TRF", period, \$4, score, period, copies
  }' assembly.trf.dat > all_repeats.bed

  awk -v pmin=${params.period_min} -v pmax=${params.period_max} -v cmin=${params.min_copies} 'BEGIN{OFS="\\t"} {
      period = \$7 + 0
      copies = \$8 + 0
      if (period >= pmin && period <= pmax && copies >= cmin) print
  }' all_repeats.bed > centromere_candidate_repeats.bed

  env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml samtools faidx ${assembly}
  cut -f1,2 ${assembly}.fai > genome.sizes

  env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml \
    bedtools makewindows -g genome.sizes -w 100000 > windows_100kb.bed
  env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml \
    bedtools coverage -a windows_100kb.bed -b centromere_candidate_repeats.bed > centromere_density_100kb.bed
  """
}

process identify_regions {
  tag "Identify highest density centromere regions"
  label 'centromere'
  publishDir "${launchDir}/results/centromere/regions", mode: 'copy', overwrite: false

  input:
  path density_bed

  output:
  path "candidate_regions.tsv", emit: regions_tsv

  script:
  """
  cut -f1 ${density_bed} | sort | uniq > scaffolds.txt

  bash ${projectDir}/../../bin/extract_candidate_region.sh ${density_bed} scaffolds.txt candidate_regions.tsv 10 ${params.region_window}
  """
}

process stringdecomposer {
  tag "${scaffold}"
  label 'centromere'
  // Closure form is required under the strict-DSL parser (Nextflow 24+): the
  // bare string was evaluated at process-definition time and failed with
  // "No such property: scaffold for class: ProcessDslV1".
  publishDir path: { "${launchDir}/results/centromere/stringdecomposer/${scaffold}" }, mode: 'copy', overwrite: false

  input:
  tuple val(scaffold), val(region_str), path(monomers)
  path assembly

  output:
  path "decomposition.tsv", emit: tsv
  path "stringdecomposer.log", emit: log
  path "${scaffold}_centromere.fa", emit: region_fa

  script:
  """
  env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml \
    samtools faidx ${assembly} "${region_str}" > ${scaffold}_centromere.fa
  if grep -q '^>' ${monomers}; then
    env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml \
      stringdecomposer ${scaffold}_centromere.fa ${monomers} -t ${task.cpus} -o decomposition_dir > stringdecomposer.log
    # the stdout above is only a log; the decomposition itself is final_decomposition.tsv
    cp decomposition_dir/final_decomposition.tsv decomposition.tsv
  else
    # no monomers (TideHunter found no tandem repeats): the aligner segfaults on an empty monomer set
    echo "skipped: no monomers" > stringdecomposer.log
    : > decomposition.tsv
  fi
  """
}

process centromere_summary {
  tag "Generate summary table"
  label 'centromere'
  publishDir "${launchDir}/results/centromere", mode: 'copy', overwrite: false

  input:
  path repeats_bed

  output:
  path "centromere_summary.tsv", emit: summary

  script:
  """
  env -u PIXI_IN_SHELL -u PIXI_PROJECT_MANIFEST -u PIXI_ENVIRONMENT_NAME -u PIXI_ENVIRONMENT_PLATFORMS -u PIXI_EXE pixi run --manifest-path ${projectDir}/../../pixi.toml \
    bedtools sort -i ${repeats_bed} > centromere_repeats_sorted.bed
  awk 'BEGIN{OFS="\\t"}
       { count[\$1]++; if(!start[\$1] || \$2<start[\$1]) start[\$1]=\$2; if(\$3>end[\$1]) end[\$1]=\$3 }
       END{ for(s in count) print s, start[s], end[s], count[s] }' \
      centromere_repeats_sorted.bed | sort -k1,1 -k2,2n > centromere_summary.tsv
  """
}

workflow CENTROMERE {
  take:
  reads
  assembly

  main:
  tidehunter_cluster(reads)
  trf_density(assembly)

  identify_regions(trf_density.out.density_bed)

  // Monomers are combined into the regions tuple (not passed as a separate
  // queue channel) and assembly arrives as a VALUE channel, so stringdecomposer
  // can fire once per region instead of being limited to a single execution.
  ch_regions = identify_regions.out.regions_tsv
    .splitCsv(sep: '\t')
    .map { row -> tuple(row[0], "${row[0]}:${row[1].toLong() + 1}-${row[2]}") }
    .combine(tidehunter_cluster.out.top_monomers)

  stringdecomposer(ch_regions, assembly)
  centromere_summary(trf_density.out.repeats_bed)
}
