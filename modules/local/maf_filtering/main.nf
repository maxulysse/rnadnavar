process MAF_FILTERING {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2:3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"

    input:
    tuple val(meta), path(maf), path(intervals)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path('*.maf'), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnadnavar/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    if [ -n "${intervals}" ]; then
        INTER=\$(awk '{print \$1 ":" \$2 "-" \$3}' ${intervals})
        interval="--interval \$INTER"
        chrom=\$(cut -f1 ${intervals} | tr "\\n" "|")
        maf=${maf}.reduced
        # This is to reduce input reading time in python in case maf is big
        grep -Ew "\${chrom}Hugo_Symbol" ${maf} > \$maf
        echo 'Reduced maf to intervals: \${chrom}'
    else
        interval=""
        maf=${maf}
    fi
    filter_mutations.py -i \$maf --output ${prefix}.maf --ref $fasta $args \$interval
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python (//;s/).*//')
    END_VERSIONS
    """

}
