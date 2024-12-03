// Define the process
process GATK_genotype_GVCFs {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "" 
	publishDir "${params.outdir}/", mode: 'copy'
  
	label 'gatk'
    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(ref_file)
	path(ref_fai)
	path(ref_dict)
	path(gvcf_file)
	path(gvcf_idx_file)


	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("${gvcf_file.simpleName}_gt.vcf"), emit: gvcf_genotyped

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	gatk --java-options "-Xmx4g" GenotypeGVCFs \
		-R ${ref_file} \
		-V ${gvcf_file}\
		-O ${gvcf_file.simpleName}_gt.vcf --include-non-variant-sites true
	"""
 }