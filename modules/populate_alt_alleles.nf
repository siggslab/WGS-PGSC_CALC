// Define the process
process populate_alt_alleles {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	publishDir "${params.outdir}/", mode: 'copy'

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(input_vcf)
	path(scoring_file)
	path(ref_file)
	path(script)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("combined_processed.vcf"), emit: combined_processed_vcf
	path("populate_alt_alleles.log"), emit: populate_alt_alleles_log

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
		#script has changed. again!
		bash $script $input_vcf "$scoring_file" $ref_file
	"""
 }
