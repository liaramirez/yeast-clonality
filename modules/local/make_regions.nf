process MAKE_REGIONS {
    label 'process_single'

    input:
    path ref

    output:
    path "regions_*.bed", emit: beds

    script:
    """
    samtools faidx ${ref}
    awk -v m=${params.min_contig_len} '\$2 >= m {print \$1"\\t0\\t"\$2}' ${ref}.fai \\
      | sort -k3,3nr > contigs.bed
    # reparte los contigs en N grupos de tamaño parecido (el mayor pendiente al grupo más liviano)
    awk -v n=${params.call_chunks} '{
          best=0; for(i=1;i<n;i++) if(load[i]<load[best]) best=i
          load[best]+=\$3; print > sprintf("regions_%03d.bed", best)
        }' contigs.bed
    """

    stub:
    """
    touch regions_000.bed regions_001.bed
    """
}
