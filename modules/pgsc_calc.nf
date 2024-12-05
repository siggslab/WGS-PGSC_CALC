// Define the process
process pgsc_calc {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	tag "" 
	publishDir "${params.outdir}/", mode: 'copy'
  	
	label 'pgsc_calc'
	container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }" 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("")

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	nextflow run /opt/pgsc_calc/main.nf \
		-profile singularity
	"""
 }