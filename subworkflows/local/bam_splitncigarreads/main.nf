//
// SPLITNCIGARREADS
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_SPLITNCIGARREADS    } from '../../../modules/nf-core/gatk4/splitncigarreads/main'
include { CRAM_MERGE_INDEX_SAMTOOLS } from '../cram_merge_index_samtools/main'

workflow BAM_SPLITNCIGARREADS {
    take:
    cram            // channel: [mandatory] [ meta, cram_markduplicates, crai ]
    dict            // channel: [mandatory] [ dict ]
    fasta           // channel: [mandatory] [ meta, fasta ]
    fasta_fai       // channel: [mandatory] [ meta, fasta_fai ]
    intervals       // channel: [mandatory] [ intervals, num_intervals ] (or [ [], 0 ] if no intervals)

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ] }

    GATK4_SPLITNCIGARREADS (
        cram_intervals,
        fasta,
        fasta_fai,
        dict.map{ meta, it ->  it  }
    )

    // Gather the recalibrated cram files
    cram_to_merge = GATK4_SPLITNCIGARREADS.out.cram.map{ meta, cram -> [ groupKey(meta, meta.num_intervals), cram ] }.groupTuple()

    // Merge and index the recalibrated cram files
    CRAM_MERGE_INDEX_SAMTOOLS(cram_to_merge, fasta, fasta_fai)

    cram_recal = CRAM_MERGE_INDEX_SAMTOOLS.out.cram_crai
    // Remove no longer necessary field: num_intervals
    .map{ meta, cram, crai -> [ meta - meta.subMap('num_intervals'), cram, crai ] }

    // Gather versions of all tools used
    versions = versions.mix(GATK4_SPLITNCIGARREADS.out.versions)
    versions = versions.mix(CRAM_MERGE_INDEX_SAMTOOLS.out.versions)


    emit:
    cram = cram_recal // channel: [ meta, cram, crai ]
    versions          // channel: [ versions.yml ]
}
