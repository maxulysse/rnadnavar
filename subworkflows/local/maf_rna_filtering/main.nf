//
// Filtering steps specific for RNA
//
include { RNA_FILTERING } from '../../../modules/local/maf_rna_filtering/main'
include { MERGE_MAF as MERGE_RNA_FILTERED_MAF } from '../../../modules/local/merge_maf/main'


workflow MAF_FILTERING_RNA {
    take:
    maf_to_filter             // maf from first pass
    maf_to_filter_realigned   // maf from realignment (second pass) [OPT?]
    fasta
    fasta_fai
    intervals
    input_sample
    realignment

    main:
    versions = Channel.empty()

    maf_to_filter.dump(tag:"rna_maf_to_filter")
    maf_to_filter_realigned.dump(tag:"rna_maf_realigned")
    final_rna_maf = Channel.empty()

    if (params.step in ['mapping', 'markduplicates', 'splitncigar',
    'prepare_recalibration', 'recalibrate', 'variant_calling', 'annotate','norm', 'consensus', 'filtering',
    'realignment', 'rna_filtering'] && (params.tools && params.tools.split(",").contains("rna_filtering"))) {

        if (params.step == 'rna_filtering') {
            maf_to_filter = input_sample
        } else {
            if (params.step == "realignment"){
                maf_to_filter = maf_to_filter_realigned // there is no first maf
                maf_to_filter_realigned = null
            }
            if (realignment){
                // RNA filtering after realignment
                maf_to_cross_first_pass = maf_to_filter
                                            .map{meta, maf -> [meta.id, meta, maf]}

                maf_to_cross_first_pass.dump(tag:"maf_to_cross_first_pass")

                maf_to_cross_realignment = maf_to_filter_realigned
                                            .map{meta, maf -> [meta.original_sample_id, meta, maf]}
                maf_to_cross_realignment.dump(tag:"maf_to_cross_realignment")
                maf_to_cross_first_pass
                                .cross(maf_to_cross_realignment).dump(tag:"maf_to_cross_crossed_pass")
                maf_crossed = maf_to_cross_first_pass
                                .cross(maf_to_cross_realignment)
                                .map { first, second ->
                                    def meta = [:]
                                    meta.patient    = first[0]
                                    meta.first_id   = first[1].id
                                    meta.second_id  = second[1].id
                                    meta.status     = first[1].status

                                    // Extract tumor_id and normal_id conditionally
                                    def id_parts_first  = first[1].id.split("_vs_")
                                    def id_parts_second = second[1].id.split("_vs_")

                                    meta.tumor_id   = id_parts_first[0]  // Always take the first part as tumor_id
                                    meta.normal_id  = id_parts_first.size() > 1 ? id_parts_first[1] : null  // Only assign normal_id if it exists

                                    // Create the new id by conditionally adding "_with_" and considering if second id contains "_vs_"
                                    def tumor_id_second = id_parts_second[0]
                                    meta.id         = meta.tumor_id + (tumor_id_second ? "_with_" + tumor_id_second : "")

                                    [meta, first[2], second[2]]
                                }

            } else {
                maf_to_filter.dump(tag:"rna_maf_to_filter2")
                maf_crossed = maf_to_filter.map{meta, maf -> [meta, maf, []]}
            }
        }
        maf_crossed.dump(tag:"maf_crossed")
        maf_crossed_intervals = maf_crossed.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for CONSENSUS module
        .map{ meta, maf, maf_realign, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], maf, maf_realign, intervals ]}

        maf_crossed_intervals.dump(tag:'maf_crossed_intervals')
        RNA_FILTERING(maf_crossed_intervals,
                    fasta,
                    fasta_fai)
        versions = versions.mix(RNA_FILTERING.out.versions)

        rna_filtering_out = RNA_FILTERING.out.maf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
        }

        rna_filtering_maf_to_merge = rna_filtering_out.intervals.map{ meta, maf ->
                                                [ groupKey(meta, meta.num_intervals), maf ]}.groupTuple()

        MERGE_RNA_FILTERED_MAF(rna_filtering_maf_to_merge)

        final_rna_maf = rna_filtering_out.no_intervals.mix(MERGE_RNA_FILTERED_MAF.out.maf_out)  //

    }


    emit:
        versions    = versions // channel: [ versions.yml ]
        maf         = final_rna_maf



}
