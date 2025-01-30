// Define the process
process index_vcf {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "INPUT: ${vcf_file}" 
	publishDir "${params.outdir}/", mode: 'copy'
	label 'gatk'
	container "${ workflow.containerEngine == 'singularity' &&
		!task.ext.singularity_pull_docker_container ?
		"${task.ext.singularity}${task.ext.singularity_version}" :
		"${task.ext.docker}${task.ext.docker_version}" }" 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(vcf_file)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("${vcf_file}.idx"), emit: vcf_idx

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
		gatk IndexFeatureFile -I ${vcf_file}
	"""
 }