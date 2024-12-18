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
	// path(bed)
	// path(bim)
	// path(fam)
	


	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("")

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
		--min_overlap 0.75 \
		--singularityCacheDir /g/data/tn36/pgs/WGS-PGSC_CALC/WGS-PGSC_CALC/singularity_cache \
		

	cp -r results/* ${workflow.projectDir}/${params.outdir}	
	"""
	//${pgsc_ref_file ? "--run_ancestry ${pgsc_ref_file} \\" : ""}
	//min_overlap=\$(cat ${min_overlap})
	//--min_overlap \$min_overlap \
	//--keep-ambiguous true

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