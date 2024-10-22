process RNA_FILTERING {
    tag "$meta.id"
    label 'process_low'

    conda null
    container 'ghcr.io/raqmanzano/rnafilt:latest'

    input:
    tuple val(meta), path(maf), path(maf_realignment), path(intervals)
    tuple val(meta1), path(fasta)
    path fasta_fai

    output:
    tuple val(meta), path('*.maf'), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnadnavar/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maf_realign_opt = maf_realignment? "true" : "false"
    """
    if [ -n "${intervals}" ]; then
        INTER=\$(awk '{print \$1 ":" \$2 "-" \$3}' ${intervals})
        interval="--interval \$INTER"
        chrom=\$(cut -f1 ${intervals} | tr "\\n" "|")
        maf=${maf}.reduced
        # This is to reduce input reading time in python in case maf is big
        grep -Ew "\${chrom}Hugo_Symbol" ${maf} > \$maf
        echo 'Reduced maf to intervals: \$chrom'
    else
        interval=""
        maf=${maf}
    fi

    if ${maf_realign_opt} ; then
        maf_realignment=${maf_realignment}.reduced
        if [ -n "${intervals}" ]; then
            # This is to reduce input reading time in python in case maf is big
            grep -Ew "\${chrom}Hugo_Symbol" ${maf_realignment} > \$maf_realignment
            echo 'Reduced realignment maf to intervals: \$chrom'
            maf_realign_opt="--maf_realign \$maf_realignment"
        else
            interval=""
        fi
    else
        maf_realign_opt=""
    fi

    filter_rna_mutations.py \\
        --maf \$maf \\
        --ref $fasta \\
        --output ${prefix}.maf \$interval \\
        \$maf_realign_opt \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python (//;s/).*//')
    END_VERSIONS
    """

}
