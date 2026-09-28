process MERGE_CALLS {
    label 'process_single'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path vcfs

    output:
    path "joint.vcf.gz",     emit: vcf
    path "joint.vcf.gz.csi", emit: csi
    path "genotypes.tsv",    emit: gt

    script:
    """
    ls regions_*.vcf.gz | sort > list.txt
    bcftools concat -f list.txt -Oz -o unsorted.vcf.gz
    bcftools sort -Oz -o joint.vcf.gz unsorted.vcf.gz
    bcftools index joint.vcf.gz
    # tabla: una fila por sitio, una columna por muestra con GT:DP
    ( bcftools query -l joint.vcf.gz | paste -sd'\\t' ;
      bcftools query -f '[%GT:%DP\\t]\\n' joint.vcf.gz | sed 's/\\t\$//' ) > genotypes.tsv
    """

    stub:
    """
    touch joint.vcf.gz joint.vcf.gz.csi
    printf "A\\tB\\tC\\n0/1:20\\t0/1:22\\t1/1:25\\n" > genotypes.tsv
    """
}
