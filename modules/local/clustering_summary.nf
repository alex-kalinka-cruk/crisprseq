process CLUSTERING_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::biopython=1.78" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/biopython:1.78' :
        'https://depot.galaxyproject.org/singularity/biopython:1.78' }"

    input:
    tuple val(meta), path(reads), path(summary)

    output:
    tuple val(meta), path(summary), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'clustering_summary.py'
}
