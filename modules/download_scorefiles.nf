process download_scorefiles {
    // labels are defined in conf/modules.config
    //label 'process_medium'
    label 'pgscatalog_utils' // controls conda, docker, + singularity options

    tag "$meta"
    time '30m'

    conda "${task.ext.conda}"

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    input:
    val(meta)
    val(build)

    output:
    path "PGS*.txt.gz" , emit: scorefiles
    path "versions.yml", emit: versions

    script:
    def accession_args = meta.pgs_id   ? "-i $meta.pgs_id"  : ""
    def traits_args = meta.efo_id   ? "-t $meta.efo_id"  : ""
    def publication_args = meta.pgp_id ? "-p $meta.pgp_id": ""

    """
    pgscatalog-download $accession_args \
        $traits_args \
        $publication_args \
        -b $build \
        -o \$PWD \
        -v \
        -c pgsc_calc/$workflow.manifest.version

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        pgscatalog.core: \$(echo \$(python -c 'import pgscatalog.core; print(pgscatalog.core.__version__)'))
    END_VERSIONS
    """
}