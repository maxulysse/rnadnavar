//
// Filtering somatic mutation analysis
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
include { MAF_FILTERING as FILTERING      } from '../../../modules/local/maf_filtering/main'
include { MERGE_MAF as MERGE_FILTERED_MAF } from '../../../modules/local/merge_maf/main'

workflow MAF_FILTERING {

    take:
    maf_to_filter
    fasta
    intervals
    input_sample
    realignment

    main:
    versions  = Channel.empty()
    maf       = Channel.empty()
    if ((params.step in ['mapping', 'markduplicates', 'splitncigar',
                'prepare_recalibration', 'recalibrate', 'variant_calling', 'annotate',
                'norm', 'consensus', 'filtering', 'realignment'] &&
                ((params.tools && params.tools.split(",").contains("filtering")))) ||
                realignment) {

        if (params.step == 'filtering') {
            maf_to_filter = input_sample
        }

        maf_to_filter.dump(tag:'maf_to_filter')
        intervals.dump(tag:'intervals_filter')

        maf_intervals = maf_to_filter.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for CONSENSUS module
        .map{ meta, maf, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], maf, intervals ]}

        maf_intervals.dump(tag:'maf_intervals_filtering')
        // BASIC FILTERING
        FILTERING(maf_intervals, fasta)

        filtering_out = FILTERING.out.maf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
        }

        filtering_maf_to_merge = filtering_out.intervals.map{ meta, maf ->
                                                [ groupKey(meta, meta.num_intervals), maf ]}.groupTuple()

        MERGE_FILTERED_MAF(filtering_maf_to_merge)

        maf      = filtering_out.no_intervals.mix(MERGE_FILTERED_MAF.out.maf_out)  // 1 consensus_maf from all callers
        versions = versions.mix(FILTERING.out.versions)
    }

    emit:
    maf        = maf
    versions   = versions                                                         // channel: [ versions.yml ]
}
