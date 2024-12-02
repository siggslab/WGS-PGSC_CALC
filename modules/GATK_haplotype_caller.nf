// Define the process
process generate_report {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "" 
	publishDir "${params.outdir}/", mode: 'copy'
  container '' 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path("")

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("")

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""

	"""
 }