// Define the process
process make_samplesheet {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "" 
	publishDir "${params.outdir}/", mode: 'copy'
  container '' 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(processed_vcf)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("samplesheet.csv"), emit: samplesheet

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	basename=\$(basename ${processed_vcf})
	cat <<EOF > samplesheet.csv
	sampleset,path_prefix,chrom,format
	cineca,\${basename},,vcf
	EOF
	"""
 }