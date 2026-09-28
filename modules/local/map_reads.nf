process MAP_READS {
    tag "$sample"
    label 'process_high'

    input:
    tuple val(sample), path(reads)
    path ref

    output:
    path "${sample}.bam",     emit: bam
    path "${sample}.bam.bai", emit: bai

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont --secondary=no \\
        -R '@RG\\tID:${sample}\\tSM:${sample}' ${ref} ${reads} \\
      | samtools sort -@ 4 -o ${sample}.bam
    samtools index ${sample}.bam
    """

    stub:
    """
    touch ${sample}.bam ${sample}.bam.bai
    """
}
