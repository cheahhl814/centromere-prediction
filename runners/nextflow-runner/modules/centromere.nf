process tidehunter_cluster {
  tag "Discover tandem repeat monomers from reads"
  label 'centromere'
  publishDir "${launchDir}/results/centromere/tidehunter", mode: 'copy', overwrite: false, pattern: 'consensus_clustered.fasta'

  input:
  path reads

  output:
  path "consensus_clustered.fasta", emit: monomers

  script:
  """
  pixi run --manifest-path ${projectDir}/../../pixi.toml \\
    TideHunter -p 50 -P 250 -c 5 -t ${task.cpus} ${reads} > consensus_raw.fasta

  pixi run --manifest-path ${projectDir}/../../pixi.toml \\
    cd-hit-est -i consensus_raw.fasta -o consensus_clustered.fasta -c 0.8 -n 5 -T ${task.cpus} -M 8000 -d 0 -sc 1 -sf 1
  """
}

process trf_density {
  tag "Map tandem repeats and compute density"
  label 'centromere'
  publishDir "${launchDir}/results/centromere/trf", mode: 'copy', overwrite: false

  input:
  path assembly

  output:
  path "centromere_candidate_repeats.bed", emit: repeats_bed
  path "centromere_density_100kb.bed", emit: density_bed
  path "genome.sizes", emit: sizes

  script:
  """
  pixi run --manifest-path ${projectDir}/../../pixi.toml \
    trf ${assembly} 2 7 7 80 10 50 250 -f -d -m -ngs -h > assembly.trf.dat 2>/dev/null || true

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

  awk -v pmin=50 -v pmax=250 'BEGIN{OFS="\\t"} {
      period = \$7 + 0
      copies = \$8 + 0
      if (period >= pmin && period <= pmax && copies >= 5) print
  }' all_repeats.bed > centromere_candidate_repeats.bed

  pixi run --manifest-path ${projectDir}/../../pixi.toml samtools faidx ${assembly}
  cut -f1,2 ${assembly}.fai > genome.sizes

  pixi run --manifest-path ${projectDir}/../../pixi.toml \
    bedtools makewindows -g genome.sizes -w 100000 > windows_100kb.bed
  pixi run --manifest-path ${projectDir}/../../pixi.toml \
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

  bash ${projectDir}/../../bin/extract_candidate_region.sh ${density_bed} scaffolds.txt candidate_regions.tsv
  """
}

process stringdecomposer {
  tag "${scaffold}"
  label 'centromere'
  publishDir "${launchDir}/results/centromere/stringdecomposer/${scaffold}", mode: 'copy', overwrite: false

  input:
  tuple val(scaffold), val(region_str), path(monomers)
  path assembly

  output:
  path "decomposition.tsv", emit: tsv
  path "${scaffold}_centromere.fa", emit: region_fa

  script:
  """
  pixi run --manifest-path ${projectDir}/../../pixi.toml \
    samtools faidx ${assembly} "${region_str}" > ${scaffold}_centromere.fa
  pixi run --manifest-path ${projectDir}/../../pixi.toml \
    stringdecomposer ${scaffold}_centromere.fa ${monomers} -t ${task.cpus} -o decomposition_dir > decomposition.tsv
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
  pixi run --manifest-path ${projectDir}/../../pixi.toml \
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
    .map { row -> tuple(row[0], "${row[0]}:${row[1]}-${row[2]}") }
    .combine(tidehunter_cluster.out.monomers)

  stringdecomposer(ch_regions, assembly)
  centromere_summary(trf_density.out.repeats_bed)
}
