#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2

// =================================================================
// main.nf is the pipeline script for a nextflow pipeline
// Should contain the following sections:
	// Process definitions
    // Channel definitions
    // Workflow structure
	// Workflow summary logs 

// Examples are included for each section. Remove them and replace
// with project-specific code. For more information see:
// https://www.nextflow.io/docs/latest/index.html.
//
// ===================================================================

// Import processes or subworkflows to be run in the workflow
// Each of these is a separate .nf script saved in modules/ directory
// See https://training.nextflow.io/basic_training/modules/#importing-modules 
include { check_input } 						from './modules/check_input'
include { generate_PRS_snp_positions_list } 	from './modules/generate_PRS_snp_positions_list.nf' 
include { remove_chr_from_BAM } 				from './modules/remove_chr_from_BAM.nf'
include { download_scorefiles } 				from './modules/download_scorefiles.nf'
include { combine_PRS_snp_positions_lists } 	from './modules/combine_PRS_snp_positions_lists.nf'
include { GATK_haplotype_caller } 				from './modules/GATK_haplotype_caller.nf'
include { GATK_genotype_GVCFs } 				from './modules/GATK_genotype_GVCFs.nf'
include { pgsc_calc } 						    from './modules/pgsc_calc.nf'
include { make_samplesheet}                     from './modules/make_samplesheet.nf'
include { calculate_min_overlap }               from './modules/calculate_min_overlap.nf'

// Print a header for your pipeline 
log.info """\

=======================================================================================
Name of the pipeline - nf 
=======================================================================================

Created by <YOUR NAME> 
Find documentation @ https://sydney-informatics-hub.github.io/Nextflow_DSL2_template_guide/
Cite this pipeline @ INSERT DOI

=======================================================================================
Workflow run parameters 
=======================================================================================
scorefile           : ${params.scorefile}
outdir             : ${params.outdir}
workDir             : ${workflow.workDir}
bamfile	            : ${params.bamfile}
target_build        : ${params.target_build}
pgs_id	        	: ${params.pgs_id}
efo_id	   	        : ${params.efo_id}
pgp_id	   	        : ${params.pgp_id}
ref 	            : ${params.ref}
dbsnp               : ${params.dbsnp}
singularityCacheDir : ${params.singularityCacheDir}
min_overlap         : ${params.min_overlap}
=======================================================================================

"""

/// Help function 
// This is an example of how to set out the help function that 
// will be run if run command is incorrect or missing. 

def helpMessage() {
    log.info"""
  Usage:  nextflow run main.nf --scorefile <scorefile.txt> -- bamfile <bamfile.bam> --outdir <output_directory>

  Required Arguments:

  --scorefile	Specify full path and name of score file.
  --bamfile		Specify full path and name of BAM file.

  Optional Arguments:

  --outdir	Specify path to output directory. 
	
""".stripIndent()
}

// Define the prepareAccessions function
def prepareAccessions(String accession, String key) {
    def unique_accession = accession.replaceAll('\\s','').tokenize(',').unique()
    def good_accessions = []
    unique_accession.each { it ->
        if (!(it ==~ "(?:PGP|PGS)[0-9]{6}|(?:[A-Za-z]+)_[0-9]+")) {
            System.err.println "WARNING: ${it} doesn't seem like a valid PGS Catalog accession, ignoring"
        } else {
            good_accessions.add(it)
        }
    }
    if (good_accessions) {
        return [(key): good_accessions.join(" ")]
    } else {
        return [(key): ""]
    }
}

// Define workflow structure. Include some input/runtime tests here.
// See https://www.nextflow.io/docs/latest/dsl2.html?highlight=workflow#workflow
workflow {

// Show help message if --help is run or (||) a required parameter (input) is not provided

if ( params.help || !params.bamfile || !params.target_build || !params.ref || !params.dbsnp ||
    (params.target_build != 'GRCh37' && params.target_build != 'GRCh38') || 
    !(params.pgs_id || params.efo_id || params.pgp_id || params.scorefile) ) {  
// Invoke the help function above and exit
	helpMessage()
	exit 1
	// consider adding some extra contigencies here.
	// could validate path of all input files in list?
	// could validate indexes for reference exist?

// If none of the above are a problem, then run the workflow
} else {
	
	// DEFINE CHANNELS 
	// See https://www.nextflow.io/docs/latest/channel.html#channels
	// See https://training.nextflow.io/basic_training/channels/ 

	// DEMO CODE: DELETE FOR YOUR OWN WORKFLOWS - VALIDATE INPUT SAMPLES 
	//check_input(Channel.fromPath(params.input, checkIfExists: true))
	
	//SUBWORKFLOW Download - Scorefiles
	ch_scores = Channel.empty()
    if (params.scorefile) {
        ch_scores = ch_scores.mix(Channel.fromPath(params.scorefile, checkIfExists: true))
    }

    // make sure accessions look sensible before querying PGS Catalog
    def pgs_id = prepareAccessions(params.pgs_id, "pgs_id")
    def pgp_id = prepareAccessions(params.pgp_id, "pgp_id")
	def efo_id = prepareAccessions(params.efo_id, "efo_id")
    
    def accessions = pgs_id + pgp_id + efo_id

    if (!accessions.every { it.value == "" }) {
        download_scorefiles(accessions, params.target_build)
        ch_scores = ch_scores.mix(download_scorefiles.out.scorefiles)
    }

    if (!params.scorefile && accessions.every { it.value == "" }) {
        Nextflow.error("No valid accessions or scoring files provided. Please double check --pgs_id, --pgp_id, --trait_efo, or --scorefile parameters")
    }

	//END SUBWORKFLOW

	//Run generate_PRS_snp_positions_list for each scorefile
	generate_PRS_snp_positions_list(ch_scores.flatten())
		.collect()
		.set { prs_snp_positions_files }
	//Combine all PRS_snp_positions.list files into one, removing duplicate lines
	combine_PRS_snp_positions_lists(prs_snp_positions_files)

    //TODO: Not sure if this is best method, should I be using GRCh37 instead?
    //This seems very resource intensive
    remove_chr_from_BAM(params.bamfile)

	//Run GATK HaplotypeCaller on the BAM file
    GATK_haplotype_caller(
        remove_chr_from_BAM.out.bam_file_no_chr,
        "${params.ref}.fasta", 
        "${params.ref}.fasta.fai",
        "${params.ref}.dict",
        combine_PRS_snp_positions_lists.out.PRS_snp_positions,
        params.dbsnp,
        "${params.dbsnp}.idx"
        )

    //Run GATK GenotypeGVCFs on the gVCF file
    GATK_genotype_GVCFs(
        "${params.ref}.fasta",
        "${params.ref}.fasta.fai",
        "${params.ref}.dict",
        GATK_haplotype_caller.out.haplotypeCalled_gvcf,
        GATK_haplotype_caller.out.haplotypeCalled_gvcf_idx,
        params.dbsnp,
        "${params.dbsnp}.idx",
        combine_PRS_snp_positions_lists.out.PRS_snp_positions
        )

    //Make samplesheet from processed VCF
    make_samplesheet(GATK_genotype_GVCFs.out.gvcf_genotyped)

    //Calculate Minimum Overlap (reference sites otherwise lower %)
    calculate_min_overlap(GATK_genotype_GVCFs.out.gvcf_genotyped, params.min_overlap)

    //Run pg_sc_calc with input files
    pgsc_calc(
        make_samplesheet.out.samplesheet,
        GATK_genotype_GVCFs.out.gvcf_genotyped,
        params.target_build,
        ch_scores.flatten().collect(),
        workflow.workDir,
        calculate_min_overlap.out.min_overlap,
        "${params.ref}.fasta", 
        "${params.ref}.fasta.fai",
        "${params.ref}.dict"
    )
}}

// Print workflow execution summary 
workflow.onComplete {
summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
results     : ${params.outdir}

=======================================================================================
  """
println summary

}
