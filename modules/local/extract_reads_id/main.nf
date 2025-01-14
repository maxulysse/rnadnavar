process SAMTOOLS_EXTRACT_READ_IDS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(input), path(index), path(bed)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*_IDs_all.txt") , emit: read_ids
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        ${reference} \\
        -L $bed \\
        $args \\
        $input | cut -f1 > ${prefix}_IDs_all.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
