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
	path(samplesheet)
	path(processed_vcf)
	val target_build
	file(scorefiles)
	path(workDir)
	val(min_overlap)
	path(singularity_cache_dir)
	


	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("results"), emit: results_dir

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	scorefile_paths=\$(printf " %s" ${scorefiles.collect { it.toString() }.join(',')})
	export SINGULARITY_CACHEDIR=${params.singularityCacheDir}
	export NXF_SINGULARITY_CACHEDIR=${params.singularityCacheDir}
	export CONDA_PKGS_DIRS=${params.singularityCacheDir}
	#I think it might be supposed to be export SINGULARITY_CACHEDIR=${params.singularityCacheDir}
	export HOME=${workDir}/home
    mkdir -p \$HOME

	export PATH=\$PATH:/opt/pbs/default/bin/:/opt/singularity/bin/singularity

	#Clear the pgsc_calc work directory as it causes weird caching issues.
	rm -rf ${workflow.projectDir}/work/pgsc_calc

	nextflow run /opt/pgsc_calc/main.nf \
		-profile conda \
		--input ${samplesheet} \
		--target_build ${target_build} \
    	--scorefile "*.txt.gz" \
		-w ${workflow.projectDir}/work/pgsc_calc \
		--outDir ${workflow.projectDir}/results \
		--min_overlap ${min_overlap} \
		--singularityCacheDir ${singularity_cache_dir} \
		

	#cp -r results/* ${workflow.projectDir}/${params.outdir}	
	"""
 }