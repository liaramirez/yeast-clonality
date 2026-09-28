process CHUNK_READS {
    tag "$sample"
    label 'process_low'

    input:
    tuple val(sample), path(fastq_dir)

    output:
    tuple val(sample), path("${sample}.reads.fq.gz"), emit: reads

    script:
    if (params.chunk_size > 0)
        """
        cat ${fastq_dir}/*.fastq.gz \\
          | seqkit sliding -W ${params.chunk_size} -s ${params.chunk_size} -g -j ${task.cpus} \\
          | gzip -1 > ${sample}.reads.fq.gz
        """
    else
        """
        cat ${fastq_dir}/*.fastq.gz > ${sample}.reads.fq.gz
        """

    stub:
    """
    echo "@r1" | gzip > ${sample}.reads.fq.gz
    """
}
