// Define the process
process generate_PRS_snp_positions_list {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "INPUT: ${scoring_file}" 
	publishDir "${params.outdir}/", mode: 'copy'
    container '' 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(scoring_file)
	path(dbsnp_file)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("${scoring_file.simpleName}_PRS_snp_positions.list"), emit: PRS_snp_positions

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	// Determine the command to read the file based on its extension
	"""
		#script has changed. again!!!!
		bash ${workflow.projectDir}/lib/generate_PRS_snp_positions_list.sh $scoring_file $dbsnp_file
	"""
 }