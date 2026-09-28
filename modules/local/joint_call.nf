process JOINT_CALL {
    tag "${bed.baseName}"
    label 'process_single'

    input:
    path bed
    path bams
    path bais
    path ref

    output:
    path "${bed.baseName}.vcf.gz", emit: vcf

    script:
    """
    samtools faidx ${ref}
    ls *.bam | sort > bams.txt
    bcftools mpileup -X ont-sup-1.20 -f ${ref} -R ${bed} -q 20 -Q 10 \\
        -a FORMAT/AD,FORMAT/DP -b bams.txt \\
      | bcftools call -mv --ploidy 2 \\
      | bcftools view -v snps -i 'QUAL>=${params.min_qual}' -Oz -o ${bed.baseName}.vcf.gz
    """

    stub:
    """
    touch ${bed.baseName}.vcf.gz
    """
}
