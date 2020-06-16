#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["atlas_directory":"$params.atlas_directory",
                "use_orientational_priors":"$params.use_orientational_priors",
                "use_bs_tracking_mask":"$params.use_bs_tracking_mask",
                "bs_tracking_mask_dilation":"$params.bs_tracking_mask_dilation",
                "use_bs_endpoints_include":"$params.use_bs_endpoints_include",
                "bs_endpoints_mask_dilation":"$params.bs_endpoints_mask_dilation",
                "use_tracking_mask_as_seeding":"$params.use_tracking_mask_as_seeding",
                "local_tracking":"$params.local_tracking",
                "pft_tracking":"$params.pft_tracking",
                "seeding":"$params.seeding",
                "nbr_seeds":"$params.nbr_seeds",
                "algo":"$params.algo",
                "basis":"$params.basis",
                "min_length":"$params.min_length",
                "max_length":"$params.max_length",
                "compress_error_tolerance":"$params.compress_error_tolerance",
                "tracking_seed":"$params.tracking_seed",
                "recobundles":"$params.recobundles",
                "wb_clustering_thr":"$params.wb_clustering_thr",
                "model_clustering_thr":"$params.model_clustering_thr",
                "prunning_thr":"$params.prunning_thr",
                "outlier_alpha":"$params.outlier_alpha"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "SCIL bundle specific pipeline"
log.info "============================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[Atlas]"
log.info "Atlas Anatomy: $params.atlas_anat"
log.info "Atlas Directory: $params.atlas_directory"
log.info ""
log.info "[Priors options]"
log.info "BS Tracking Mask: $params.use_bs_tracking_mask"
log.info "BS Endpoints Mask: $params.use_bs_endpoints_include"
log.info "Endpoints Mask Dilation: $params.bs_endpoints_mask_dilation"
log.info "Seeding From Tracking Mask: $params.use_tracking_mask_as_seeding"
log.info "Tracking Mask Dilation: $params.bs_tracking_mask_dilation"
log.info ""
log.info "[Tracking options]"
log.info "Local Tracking: $params.local_tracking"
log.info "PFT Tracking: $params.pft_tracking"
log.info "Algo: $params.algo"
log.info "Seeding Type: $params.seeding"
log.info "Number of Seeds: $params.nbr_seeds"
log.info "Random Seed: $params.tracking_seed"
log.info "Minimum Length: $params.min_length"
log.info "Maximum Length: $params.max_length"
log.info "Compressing Threshold: $params.compress_error_tolerance"
log.info ""
log.info "[Recobundles options]"
log.info "Segmentation with Recobundles: $params.recobundles"
log.info "Whole Brain Clustering Threshold: $params.wb_clustering_thr"
log.info "Model Clustering Threshold: $params.model_clustering_thr"
log.info "Prunning Threshold: $params.prunning_thr"
log.info "Outlier Removal Alpha: $params.outlier_alpha"
log.info ""
log.info ""

log.info "Input: $params.root"
root = file(params.root)
/* Watch out, files are ordered alphabetically in channel */
    in_data = Channel
        .fromFilePairs("$root/**/{*fa.nii.gz,*fodf.nii.gz,*tracking_mask.nii.gz}",
                        size: 3,
                        maxDepth:2,
                        flat: true) {it.parent.name}

    map_pft = Channel
        .fromFilePairs("$root/**/{*map_exclude.nii.gz,*map_include.nii.gz}",
                        size: 2,
                        maxDepth:2,
                        flat: true) {it.parent.name}

    atlas_anat = Channel.fromPath("$params.atlas_anat")
    atlas_bundles = Channel.fromPath("$params.atlas_directory/*.trk")
    algo_list = params.algo?.tokenize(',')

(anat_for_registration, anat_for_deformation, fod_and_mask_for_priors) = in_data
    .map{sid, anat, fodf, tracking_mask -> 
        [tuple(sid, anat),
        tuple(sid, anat),
        tuple(sid, fodf, tracking_mask)]}
    .separate(3)

(map_in_for_tracking, map_ex_for_tracking) = map_pft
    .map{sid, map_exclude , map_include -> 
        [tuple(sid, map_include),
        tuple(sid, map_exclude)]}
    .separate(2)

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

anat_for_registration
    .combine(atlas_anat)
    .set{anats_for_registration}
process Register_Anat {
    cpus params.register_processes
    input:
    set sid, file(native_anat), file(atlas) from anats_for_registration

    output:
    set sid, "${sid}__output1InverseWarp.nii.gz", "${sid}__output0GenericAffine.mat" into deformation_for_warping
    file "${sid}__outputWarped.nii.gz"
    script:
    """
    antsRegistrationSyNQuick.sh -d 3 -f ${native_anat} -m ${atlas} -n ${params.register_processes} -o ${sid}__output
    """ 
}


anat_for_deformation
    .join(deformation_for_warping)
    .set{anat_deformation_for_warp}
process Warp_Bundle {
    cpus 2
    input:
    set sid, file(anat), file(warp), file(affine) from anat_deformation_for_warp
    each file(bundle_name) from atlas_bundles

    output:
    set sid, val(bundle_name.baseName), "${sid}__${bundle_name.baseName}_warp.trk" into bundles_for_priors, models_for_recobundles
    script:
    """
    ConvertTransformFile 3 ${affine} ${affine}.txt --hm --ras
    scil_apply_transform_to_tractogram.py ${bundle_name} ${warp} ${affine}.txt ${bundle_name.baseName}_linear.trk --inverse --remove_invalid -f
    scil_apply_warp_to_tractogram.py ${bundle_name.baseName}_linear.trk ${anat} ${warp} ${sid}__${bundle_name.baseName}_warp.trk --remove_invalid -f
    """ 
}


fod_and_mask_for_priors
    .combine(bundles_for_priors, by: 0)
    .set{fod_mask_bundles_for_priors}

process Generate_Priors {
    cpus 2
    errorStrategy 'ignore'
    publishDir = {"./results_bst/$sid/$task.process/${bundle_name}"}
    input:
    set sid, file(fod), file(mask), val(bundle_name), file(bundle) from fod_mask_bundles_for_priors

    output:
    set sid, val(bundle_name), "${sid}__${bundle_name}_efod.nii.gz" into efod_for_tracking
    set sid, val(bundle_name), "${fod}" into fod_for_tracking
    set sid, "${sid}__${bundle_name}_priors.nii.gz"
    set sid, val(bundle_name), "${sid}__${bundle_name}_todi_mask_dilate.nii.gz", \
        "${sid}__${bundle_name}_endpoints_mask_dilate.nii.gz" into masks_for_seeding
    set sid, val(bundle_name), "${mask}", "${sid}__${bundle_name}_todi_mask_dilate.nii.gz" into masks_for_tracking
    set sid, val(bundle_name), "${sid}__${bundle_name}_todi_mask_dilate.nii.gz" into masks_for_map_ex
    set sid, val(bundle_name), "${sid}__${bundle_name}_endpoints_mask_dilate.nii.gz" into masks_for_map_in
    script:
    """
    scil_generate_priors_from_bundle.py ${bundle} ${fod} ${mask} \
        --sh_basis $params.basis --output_prefix ${sid}__${bundle_name}_
    scil_image_math.py dilation ${sid}__${bundle_name}_todi_mask.nii.gz \
        $params.bs_tracking_mask_dilation dilate_todi.nii.gz --data_type uint8
    scil_image_math.py multiplication ${mask} dilate_todi.nii.gz \
        ${sid}__${bundle_name}_todi_mask_dilate.nii.gz --data_type uint8

    scil_image_math.py dilation ${sid}__${bundle_name}_endpoints_mask.nii.gz \
        $params.bs_endpoints_mask_dilation dilate_endpoints.nii.gz --data_type uint8
    scil_image_math.py multiplication ${mask} dilate_endpoints.nii.gz \
        ${sid}__${bundle_name}_endpoints_mask_dilate.nii.gz --data_type uint8
    """
}

process Seeding_Mask {
    cpus 2
    input:
    set sid, val(bundle_name), file(tracking_mask), file(endpoints_mask) from masks_for_seeding

    output:
    set sid, val(bundle_name), "${sid}__${bundle_name}_seeding_mask.nii.gz" into \
        seeding_mask_for_PFT_tracking, seeding_mask_for_local_tracking
    script:
    if (params.use_tracking_mask_as_seeding)
        """
        mv ${tracking_mask} ${sid}__${bundle_name}_seeding_mask.nii.gz
        """
    else
        """
        scil_image_math.py multiplication ${tracking_mask} ${endpoints_mask} \
            ${sid}__${bundle_name}_seeding_mask.nii.gz --data_type uint8
        """
}

process Tracking_Mask {
    cpus 2
    input:
    set sid, val(bundle_name), file(tracking_mask), file(bs_mask) from masks_for_tracking

    output:
    set sid, val(bundle_name), "${sid}__${bundle_name}_tracking_mask.nii.gz" \
        into tracking_mask_for_local_tracking
    when: 
    params.local_tracking
    script:
    if (params.use_bs_tracking_mask)
        """
        mv ${bs_mask} ${sid}__${bundle_name}_tracking_mask.nii.gz
        """
    else
        """
        mv ${tracking_mask} ${sid}__${bundle_name}_tracking_mask.nii.gz
        """
}

if (params.use_orientational_priors)
    efod_for_tracking
        .into{fod_for_local_tracking; fod_for_PFT_tracking} 
else
    fod_for_tracking
        .into{fod_for_local_tracking; fod_for_PFT_tracking}
tracking_mask_for_local_tracking
    .combine(fod_for_local_tracking, by: [0,1])
    .combine(seeding_mask_for_local_tracking, by: [0,1])
    .set{mask_seeding_mask_fod_for_tracking}
process Local_Tracking {
    cpus 2
    input:
    set sid, val(bundle_name), file(tracking_mask), file(efod), file(seeding_mask) \
        from mask_seeding_mask_fod_for_tracking
    each algo from algo_list

    output:
    set sid, val(bundle_name), val(algo), val('local'), \
        "${sid}__${bundle_name}_${algo}_${params.seeding}_${params.nbr_seeds}.trk" into \
            local_bundles_for_recobundles
    when: 
    params.local_tracking
    script:
    """
    scil_compute_local_tracking.py ${efod} ${seeding_mask} ${tracking_mask} \
        ${sid}__${bundle_name}_${algo}_${params.seeding}_${params.nbr_seeds}.trk \
        --sh_basis $params.basis --min_len $params.min_length --max_len $params.max_length \
        --$params.seeding $params.nbr_seeds --compress $params.compress_error_tolerance \
        --seed $params.tracking_seed --algo ${algo}
    """
}

map_in_for_tracking
    .combine(masks_for_map_in, by: 0)
    .set{masks_map_in_for_bs}
process Generate_Map_Include {
    cpus 2
    input:
    set sid, file(map_include), val(bundle_name), file(endpoints_mask) from masks_map_in_for_bs

    output:
    set sid, val(bundle_name), "${sid}__${bundle_name}_map_include.nii.gz" into map_in_for_PFT_tracking
    when: 
    params.pft_tracking
    script:
    if (params.use_bs_endpoints_include)
        """
        scil_image_math.py dilation ${endpoints_mask} $params.bs_endpoints_mask_dilation \
            dilate_endpoints.nii.gz --data_type uint8
        scil_image_math.py multiplication dilate_endpoints.nii.gz ${map_include} \
            ${sid}__${bundle_name}_map_include.nii.gz --data_type float32
        """
    else
        """
        mv $map_include ${sid}__${bundle_name}_map_include.nii.gz
        """ 
}

map_ex_for_tracking
    .combine(masks_for_map_ex, by: 0)
    .set{masks_map_ex_for_bs}
process Generate_Map_Exclude {
    cpus 2
    input:
    set sid, file(map_exclude), val(bundle_name), file(tracking_mask) from masks_map_ex_for_bs

    output:
    set sid, val(bundle_name), "${sid}__${bundle_name}_map_exclude.nii.gz" into \
        map_ex_for_PFT_tracking
    when: 
    params.pft_tracking
    script:
    if (params.use_bs_tracking_mask)
        """
        scil_image_math.py invert ${tracking_mask} inverted_mask.nii.gz \
            --data_type uint8
        
        scil_image_math.py addition ${map_exclude} inverted_mask.nii.gz \
            ${sid}__${bundle_name}_map_exclude.nii.gz --data_type float32
        """
    else
        """
        mv ${map_exclude} ${sid}__${bundle_name}_map_exclude.nii.gz
        """
}

map_ex_for_PFT_tracking
    .combine(map_in_for_PFT_tracking, by: [0,1])
    .combine(fod_for_PFT_tracking, by: [0,1])
    .combine(seeding_mask_for_PFT_tracking, by: [0,1])
    .set{maps_seeding_mask_fod_for_tracking}
process PFT_Tracking {
    cpus 2
    input:
    set sid, val(bundle_name), file(map_exclude), file(map_include), file(efod), \
        file(seeding_mask) from maps_seeding_mask_fod_for_tracking
    each algo from algo_list

    output:
    set sid, val(bundle_name), val(algo), val('pft'), \
        "${sid}__${bundle_name}_${algo}_${params.seeding}_${params.nbr_seeds}.trk" into \
        pft_bundles_for_recobundles
    when: 
    params.pft_tracking
    script:
    seeding = params.seeding == 'nts' ? 'nt' : params.seeding
    """
    scil_compute_pft.py ${efod} ${seeding_mask} ${map_include} ${map_exclude} \
        ${sid}__${bundle_name}_${algo}_${params.seeding}_${params.nbr_seeds}.trk \
        --algo $algo --sh_basis $params.basis --min_length $params.min_length \
        --max_length $params.max_length --$seeding $params.nbr_seeds \
        --compress $params.compress_error_tolerance --seed $params.tracking_seed
    scil_remove_invalid_streamlines.py ${sid}__${bundle_name}_${algo}_${params.seeding}_${params.nbr_seeds}.trk  ${sid}__${bundle_name}_${algo}_${params.seeding}_${params.nbr_seeds}.trk -f
    """
}

local_bundles_for_recobundles
    .concat(pft_bundles_for_recobundles)
    .set{bundles_for_recobundles}

bundles_for_recobundles
    .combine(models_for_recobundles, by: [0,1])
    .set{bundles_models_for_recobundles}
process Recobundles_Segmentation {
    cpus 2
    publishDir = {"./results_bst/$sid/$task.process/"}
    input:
    set sid, val(bundle_name), val(algo), val(tracking_source), file(bundle), file(model) from \
        bundles_models_for_recobundles

    output:
    set sid, val(bundle_name), val(algo), val(tracking_source), "${sid}__${bundle_name}_${algo}_${tracking_source}_segmented.trk" into bundles_for_outliers
    when: 
    params.recobundles
    script:
    """
    printf "1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1" >> identity.txt
    scil_remove_invalid_streamlines.py ${bundle} tmp.trk --remove_single_point
    scil_recognize_single_bundle.py tmp.trk ${model} identity.txt \
        ${sid}__${bundle_name}_${algo}_${tracking_source}_segmented.trk \
        --tractogram_clustering_thr $params.wb_clustering_thr \
        --model_clustering_thr $params.model_clustering_thr \
        --slr_threads 1 --pruning_thr $params.prunning_thr
    """
}

process Outliers_Removal {
    cpus 2
    publishDir = {"./results_bst/$sid/$task.process/"}
    errorStrategy 'ignore'
    input:
    set sid, val(bundle_name), val(algo), val(tracking_source), file(bundle) from \
        bundles_for_outliers

    output:
    file "${sid}__${bundle_name}_${algo}_${tracking_source}_cleaned.trk"
    when: 
    params.recobundles
    script:
    """
    scil_detect_streamlines_loops.py ${bundle} no_loops.trk -a 300
    scil_outlier_rejection.py no_loops.trk  \
        ${sid}__${bundle_name}_${algo}_${tracking_source}_cleaned.trk \
        --alpha $params.outlier_alpha
    """
}
