// Define the process
process calculate_min_overlap {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	// debug = false //turn to true to print command stdout to screen
	tag "" 
	publishDir "${params.outdir}/", mode: 'copy'
	container '' 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(genotyped_vcf)
	val min_overlap_param

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("min_overlap.txt"), emit: min_overlap

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	# Bash script to calculate the minimum overlap
	# Count lines where the ALT column is empty (indicating a reference site)
	
	
	reference_sites_count=\$(grep -v "^#" ${genotyped_vcf} | awk '\$5 == "." {count++} END {print count}')
	total_sites_count=\$(grep -v "^#" ${genotyped_vcf} | wc -l)
	
	min_overlap=\$(echo "scale=10; ((\${total_sites_count}-\${reference_sites_count})/\${total_sites_count}*${min_overlap_param})" | bc)
	echo \$min_overlap > min_overlap.txt
	"""
 }