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
	//container "pgsc_calc_amd64.sif"

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(samplesheet)
	path(processed_vcf)
	val target_build
	file(scorefiles)
	path(workDir)
	path(min_overlap)
	


	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("")

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	scorefile_paths=\$(printf " %s" ${scorefiles.collect { it.toString() }.join(',')})
	export NXF_SINGULARITY_CACHEDIR=${params.singularityCacheDir}
	min_overlap=\$(cat ${min_overlap})
	
	export HOME=${workDir}/home
    mkdir -p \$HOME

	nextflow run /opt/pgsc_calc/main.nf \
		-profile conda \
		--input ${samplesheet} \
		--target_build ${target_build} \
    	--scorefile "*.txt.gz" \
		-w ${workflow.projectDir}/work/pgsc_calc \
		--min_overlap \$min_overlap \
		--outDir ${workflow.projectDir}/results \

	cp -r results/* ${workflow.projectDir}/${params.outdir}	
	"""
	// export QUARTO_CACHE=${workDir}/quarto_cache
    // mkdir -p \$QUARTO_CACHE

	//\${scorefile_paths}
	/*"""
	cd ${workflow.projectDir}/work/pgsc_calc
	echo "Workflow project directory: ${workflow.projectDir}"
    echo "Work directory: ${workflow.projectDir}/work/pgsc_calc"
	nextflow run ${workflow.projectDir}/pgsc_calc_workflow/main.nf \
	-profile test,singularity \
	--workDir ${workflow.projectDir}/work/pgsc_calc \
	-with-trace ${workflow.projectDir}/work/pgsc_calc
	"""*/
	/*"""
	scorefile_paths=\$(printf " %s" ${scorefiles.collect { it.toString() }.join(' ')})
    nextflow run /opt/pgsc_calc/main.nf \
        -profile singularity \
        --input ${samplesheet} \
        --target_build ${target_build} \
        --scorefile \${scorefile_paths}
	"""*/
 }