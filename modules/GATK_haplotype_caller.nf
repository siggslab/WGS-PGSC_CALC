// Define the process
process GATK_haplotype_caller {
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
	path(bam_file)
	path(ref_file)
	path(ref_fai)
	path(ref_dict)
	path(prs_snp_pos_file)
	path(dbsnp_file)
	path(dbsnp_idx_file)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("${bam_file.simpleName}-haplotypeCalled.g.vcf"), emit: haplotypeCalled_gvcf
	path("${bam_file.simpleName}-haplotypeCalled.g.vcf.idx"), emit: haplotypeCalled_gvcf_idx

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	gatk --java-options "-Xmx4G" HaplotypeCaller \
		-R ${ref_file} \
		-L ${prs_snp_pos_file} \
		-I ${bam_file}  \
		-O ${bam_file.simpleName}-haplotypeCalled.g.vcf -ERC GVCF \
		--dbsnp ${dbsnp_file}
	"""
 }