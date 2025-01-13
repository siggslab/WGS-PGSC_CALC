// Define the process
process create_ref_idx {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "INPUT: ${ref_file}" 
	publishDir "${params.outdir}/", mode: 'copy'
	label 'samtools'
	container "${ workflow.containerEngine == 'singularity' &&
		!task.ext.singularity_pull_docker_container ?
		"${task.ext.singularity}${task.ext.singularity_version}" :
		"${task.ext.docker}${task.ext.docker_version}" }" 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(ref_file)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("${ref_file.simpleName}.fasta.fai"), emit: ref_idx

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
		samtools faidx ${ref_file}
	"""
 }