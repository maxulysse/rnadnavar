process MAF2BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path('*.bed') , emit: bed
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env bash
    awk 'BEGIN {
        FS=OFS="\t";
        header=1;
    }
    NR==1 {
        for (i=1; i<=NF; i++) {
            if (\$i == "Chromosome") chrom_col = i;
            if (\$i == "Start_Position") start_col = i;
            if (\$i == "End_Position") end_col = i;
        }
        next;
    }
    !/^#/ {
        print \$chrom_col, \$start_col, \$end_col;
    }' ${maf} > ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$( awk --version 2>&1 | head -n 1 )
    END_VERSIONS
    """
}
