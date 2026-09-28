#!/usr/bin/env nextflow
/*
 * yeast-clonality
 * Genotype-based clonality assessment for highly heterozygous yeast genomes
 * sequenced with Oxford Nanopore.
 *
 * Steps:
 *   1. CHUNK_READS    split reads to a fixed length (read-length normalization)
 *   2. MAP_READS      map every sample to one common reference
 *   3. MAKE_REGIONS   keep contigs >= min_contig_len and split them into chunks
 *   4. JOINT_CALL     joint genotyping of all samples, one job per region chunk
 *   5. MERGE_CALLS    concatenate the per-chunk VCFs
 *   6. CLONALITY      pairwise discordance, opposite homozygotes, group splits
 */

nextflow.enable.dsl = 2

include { CHUNK_READS  } from './modules/local/chunk_reads'
include { MAP_READS    } from './modules/local/map_reads'
include { MAKE_REGIONS } from './modules/local/make_regions'
include { JOINT_CALL   } from './modules/local/joint_call'
include { MERGE_CALLS  } from './modules/local/merge_calls'
include { CLONALITY    } from './modules/local/clonality'

def helpMessage() {
    log.info """
    yeast-clonality

    Usage:
      nextflow run main.nf -profile leftraru \\
          --samplesheet samplesheet.csv \\
          --reference   ref.fa \\
          --outdir      results

    Required:
      --samplesheet   CSV with columns: sample,fastq_dir
      --reference     FASTA of the common reference (one of the samples' assemblies)

    Optional:
      --chunk_size      Split reads into pieces of this length (bp); 0 = no split [${params.chunk_size}]
      --min_contig_len  Only genotype contigs at least this long (bp)              [${params.min_contig_len}]
      --call_chunks     Number of parallel genotyping jobs                          [${params.call_chunks}]
      --min_qual        Minimum SNP QUAL                                            [${params.min_qual}]
      --min_dp          Minimum per-sample depth for a genotype                     [${params.min_dp}]
      --max_dp_factor   Maximum depth as a multiple of the sample's median depth    [${params.max_dp_factor}]
      --group           Comma-separated samples to test for sub-groups (e.g. a clone) [none]
      --outdir          Output directory                                            [${params.outdir}]
    """.stripIndent()
}

workflow {

    if (params.help) { helpMessage(); exit 0 }
    if (!params.samplesheet) error "Missing --samplesheet"
    if (!params.reference)   error "Missing --reference"

    ref = file(params.reference, checkIfExists: true)

    samples = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.fastq_dir) error "Samplesheet needs columns: sample,fastq_dir"
            tuple(row.sample, file(row.fastq_dir, checkIfExists: true))
        }

    CHUNK_READS(samples)
    MAP_READS(CHUNK_READS.out.reads, ref)

    MAKE_REGIONS(ref)
    regions = MAKE_REGIONS.out.beds.flatten()

    bams = MAP_READS.out.bam.collect()
    bais = MAP_READS.out.bai.collect()
    JOINT_CALL(regions, bams, bais, ref)

    MERGE_CALLS(JOINT_CALL.out.vcf.collect())
    CLONALITY(MERGE_CALLS.out.gt)
}

workflow.onComplete {
    log.info(workflow.success
        ? "\nDone. Results in: ${params.outdir}\n"
        : "\nFailed. Check .nextflow.log and the work/ directory of the failed task.\n")
}
