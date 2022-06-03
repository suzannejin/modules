process TIDDIT_SV {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::tiddit=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiddit:3.0.0--py39h59fae87_1' :
        'quay.io/biocontainers/tiddit:3.0.0--py39h59fae87_1' }"

    input:
    tuple val(meta), path(input)
    path  fasta
    path  fai
    path  bwa_index

    output:
    tuple val(meta), path("*.vcf")        , emit: vcf
    tuple val(meta), path("*.ploidy.tab") , emit: ploidy
    tuple val(meta), path("*.signals.tab"), emit: signals
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--ref $fasta" : ""
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    tiddit \\
        --sv \\
        $args \\
        --bam $input \\
        $reference \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.ploidy.tab
    touch ${prefix}.signals.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
    END_VERSIONS
    """
}
