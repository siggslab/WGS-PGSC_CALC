// Define the process
process remove_chr_from_BAM {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "INPUT: ${bam_file}" 
	publishDir "${params.outdir}/", mode: 'copy'

    label 'samtools'
    //container 'quay.io/biocontainers/samtools:1.14--hb421002_0' 
    //TODO change this to dynamic container
    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
    path(bam_file)

    output:
    path("${bam_file.baseName}_no_chr.bam"), emit: bam_file_no_chr

    script:
    """
    base_name=\$(basename ${bam_file} .bam)
    samtools view -h ${bam_file} | sed 's/chr//g' | samtools view -Shb - -o \${base_name}_no_chr.bam
    samtools index \${base_name}_no_chr.bam
    """
 }