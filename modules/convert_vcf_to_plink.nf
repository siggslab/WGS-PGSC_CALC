// Define the process
process convert_vcf_to_plink {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "" 
	publishDir "${params.outdir}/", mode: 'copy'
  	
	label 'plink2'
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
	path("${vcf_file.baseName}.bed"), emit: bed
	path("${vcf_file.baseName}.bim"), emit: bim
	path("${vcf_file.baseName}.fam"), emit: fam

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
		plink2 --vcf ${vcf_file} --make-bed --out ${vcf_file.baseName}
	"""
 }