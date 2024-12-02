// Define the process
process remove_chr_from_BAM {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "INPUT: ${bam_file}" 
	publishDir "${params.outdir}/", mode: 'copy'
    container 'biocontainers/samtools' 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(bam_file)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("bam_file_no_chr.bam"), emit: bam_file_no_chr

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	samtools view -h /g/data/tn36/pgs/data/wgs-low-cov/GFMC0034/GFMC0034_sorted.bam | sed 's/chr//g' | samtools view -Shb - -o GFMC0034_sorted_nochr.bam
	samtools index GFMC0034_sorted_nochr.bam
	"""
 }