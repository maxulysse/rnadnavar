//
// TUMOR ONLY VARIANT CALLING
// Should be only run on patients without normal sample
//

include { BAM_VARIANT_CALLING_SINGLE_STRELKA          } from '../bam_variant_calling_single_strelka/main'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_SAGE         } from '../bam_variant_calling_tumor_only_sage/main'
include { BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2      } from '../bam_variant_calling_tumor_only_mutect2/main'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_ALL {
    take:
    tools                         // Mandatory, list of tools to apply
    cram                          // channel: [mandatory] cram
    dict                          // channel: [mandatory] dict
    fasta                         // channel: [mandatory] fasta
    fasta_fai                     // channel: [mandatory] fasta_fai
    germline_resource             // channel: [optional]  germline_resource
    germline_resource_tbi         // channel: [optional]  germline_resource_tbi
    intervals                     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    intervals_bed_gz_tbi          // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
    intervals_bed_gz_tbi_combined // channel: [mandatory] intervals/target regions in one file zipped
    panel_of_normals              // channel: [optional]  panel_of_normals
    panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi
    joint_mutect2                 // boolean: [mandatory] [default: false] run mutect2 in joint
    realignment                   // boolean: [mandatory] [default: false] force tools to run due to realignment mode
    no_intervals

    main:
    versions = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    vcf_mutect2     = Channel.empty()
    vcf_strelka     = Channel.empty()
    vcf_sage        = Channel.empty()

    contamination_table_mutect2 = Channel.empty()
    segmentation_table_mutect2  = Channel.empty()
    artifact_priors_mutect2     = Channel.empty()

    cram.dump(tag:'tumour_only_cram')
    fasta.dump(tag:'fasta_HERE')
    // SAGE
    if (tools && tools.split(',').contains('sage') || realignment) {

        BAM_VARIANT_CALLING_TUMOR_ONLY_SAGE(
            cram,
            dict,
            fasta,
            fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] },
            intervals)

        vcf_sage = BAM_VARIANT_CALLING_TUMOR_ONLY_SAGE.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_SAGE.out.versions)
    }

    // STRELKA
    if (tools && tools.split(',').contains('strelka') || realignment) {

        BAM_VARIANT_CALLING_SINGLE_STRELKA(
            cram,
            dict,
            fasta.map{ meta, fasta ->  fasta  },
            fasta_fai,
            intervals_bed_gz_tbi,
            no_intervals
        )

        vcf_strelka = BAM_VARIANT_CALLING_SINGLE_STRELKA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SINGLE_STRELKA.out.versions)
    }

    // MUTECT2
    if (tools && tools.split(',').contains('mutect2') || realignment) {
        cram_mutect = cram.map{ meta, cram, crai ->
                joint_mutect2 ?
                //we need to keep all fields and then remove on a per-tool-basis to ensure proper joining at the filtering step
                [ meta - meta.subMap('data_type') + [ id:meta.patient ], cram, crai ] :
                [ meta - meta.subMap('data_type'), cram, crai ]
            }
        BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2(
            // Adjust meta.map to simplify joining channels
            cram_mutect,
            fasta,
            fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] },
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals,
            joint_mutect2,
            realignment
        )

        vcf_mutect2                 = BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.vcf_filtered
        contamination_table_mutect2 = BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.contamination_table
        segmentation_table_mutect2  = BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.segmentation_table
        artifact_priors_mutect2     = BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.artifact_priors
        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2.out.versions)
    }

    vcf_mutect2.dump(tag:'vcf_mutect')
    vcf_mutect2.dump(tag:'vcf_strelka')
    vcf_mutect2.dump(tag:'vcf_sage')

    vcf_all = Channel.empty().mix(
        vcf_mutect2,
        vcf_strelka,
        vcf_sage
    )

    emit:
    vcf_all
    vcf_sage
    vcf_mutect2
    vcf_strelka
    contamination_table_mutect2 = contamination_table_mutect2
    segmentation_table_mutect2  = segmentation_table_mutect2
    artifact_priors_mutect2     = artifact_priors_mutect2

    versions = versions
}
