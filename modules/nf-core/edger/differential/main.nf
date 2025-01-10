process EDGER_DIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ad/ad351ac3cb501101858fa44f1ac097f292c61228216169c39eeb249b7ed29422/data':
        'community.wave.seqera.io/library/bioconductor-edger:4.4.0--fb6573efa683e367' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(counts)

    output:
    tuple val(meta), path("*.edger.results.tsv"), emit: results
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'edger_differential.R'

    stub:
    """
    touch ${meta.id}.edger.results.tsv
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """
}
