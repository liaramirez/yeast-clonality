# yeast-clonality

Nextflow pipeline to assess clonality and strain relatedness in highly heterozygous diploid yeasts from Oxford Nanopore data.

## Why

In highly heterozygous yeasts (e.g. *Metschnikowia pulcherrima*, ~3% heterozygosity), common shortcuts mislead:

- **Assembly ANI overstates divergence.** Each assembly picks haplotypes at random along the genome, so two assemblies of the *same* clone can differ by ~1.5%.
- **Read length creates false genotype differences.** Samples with longer reads resolve heterozygous sites that short-read samples call homozygous, producing spurious "sub-groups".

This pipeline compares **genotypes**, not assemblies: all samples are genotyped jointly on one common reference after normalizing read length, and pairs are scored by **opposite homozygotes** (0/0 vs 1/1), which haplotype mosaicism and read-length effects cannot create.

## Steps

| Step | What it does |
| --- | --- |
| `CHUNK_READS` | Splits reads into fixed-length pieces (`--chunk_size`, default 2 kb) |
| `MAP_READS` | Maps every sample to the common reference (minimap2) |
| `MAKE_REGIONS` | Keeps contigs ≥ `--min_contig_len` and splits them into `--call_chunks` groups |
| `JOINT_CALL` | Joint genotyping of all samples per region group (bcftools, ONT sup preset) |
| `MERGE_CALLS` | Concatenates per-region VCFs and writes a genotype table |
| `CLONALITY` | Pairwise discordance, opposite homozygotes and, for `--group`, 2-vs-rest splits |

## Requirements

Nextflow ≥ 23.04 and, on the path of each task: `minimap2`, `samtools`, `bcftools ≥ 1.20`, `seqkit`, `python3` (see `environment.yml`).

## Quick start

1. Make a samplesheet (see `assets/samplesheet_example.csv`):

   ```csv
   sample,fastq_dir
   Lev_1,/path/to/fastq_pass/barcode01
   Lev_2,/path/to/fastq_pass/barcode02
   ```

   Each `fastq_dir` must contain that sample's `*.fastq.gz` files.

2. Run on a SLURM cluster using an existing micromamba environment:

   ```bash
   nextflow run main.nf -profile leftraru \
       --samplesheet samplesheet.csv \
       --reference   Lev_6.fa \
       --group       'Lev_1,Lev_2,Lev_3,Lev_4,Lev_5' \
       --conda_env   ani \
       --outdir      results
   ```

3. Structural test (no tools needed, seconds):

   ```bash
   nextflow run main.nf -profile test -stub
   ```

## Outputs

```
results/
├── variants/   joint.vcf.gz, genotypes.tsv
├── clonality/  pairwise.tsv, group_splits.tsv, summary.txt
└── pipeline_info/  report.html, timeline.html
```

**Reading `pairwise.tsv`:** clonal isolates show ~0 opposite homozygotes and discordance at the noise level; distinct strains show thousands of opposite homozygotes.

**Reading `group_splits.tsv`:** a real sub-lineage appears as one split well above the others; if all splits are similar, there is no sub-structure.

## Choosing the reference

Use the assembly of a sample with good coverage that is **not** unusually fragmented or haplotype-duplicated. Results should not depend on this choice; rerunning with a second reference is a useful check.

## Status

Version 0.1.0 — core clonality workflow. Planned: rDNA variant quantification, mitochondrial comparison, read simulation for benchmarking, containers.
