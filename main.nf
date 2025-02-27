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
include { convert_vcf_to_plink }                from './modules/convert_vcf_to_plink.nf'
include { bcftools_normalise_vcf }              from './modules/bcftools_normalise_vcf.nf'
include {correct_alt_for_homref}                from './modules/correct_alt_for_homref.nf'
include {populate_alt_alleles}                  from './modules/populate_alt_alleles.nf'
include {index_vcf}                             from './modules/index_vcf.nf'
include {create_ref_dict}                       from './modules/create_ref_dict.nf'
include {create_ref_idx}                        from './modules/create_ref_idx.nf'

// Print a header for your pipeline 
log.info """\

=======================================================================================
WGS PGSC_CALC - nf 
=======================================================================================

Created by Elijah Bradford 

=======================================================================================
Workflow run parameters 
=======================================================================================
scorefile           : ${params.scorefile}
outdir              : ${params.outdir}
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
run_ancestry        : ${params.run_ancestry}
add_sex             : ${params.add_sex}
=======================================================================================

"""

/// Help function 

def helpMessage() {
    log.info"""
  Usage:  nextflow run main.nf --bamfile <bamfile.bam> --target_build <GRCh37|GRCh38> --ref <reference.fasta> --dbsnp <dbsnp.vcf> [--scorefile <scorefile.txt> | --pgs_id <pgs_id> | --efo_id <efo_id> | --pgp_id <pgp_id>] --outdir <output_directory>

  Required Arguments:

  --bamfile                 Specify full path and name of BAM file.
  --target_build            Specify the target build. Must be 'GRCh37' or 'GRCh38'.
  --ref                     Specify full path and name of the reference genome file.
  --dbsnp                   Specify full path and name of the dbSNP file.
  --outdir                  Specify path to output directory.
  
  At least one of the following must be provided:
  --scorefile               Specify full path and name of score file.
  --pgs_id                  Specify PGS Catalog ID.
  --efo_id                  Specify EFO ID.
  --pgp_id                  Specify PGP ID.

  Optional Arguments:

  --help                    Show this help message and exit.
  --singularityCacheDir     Specify path to Singularity cache directory.
  --min_overlap             Specify minimum overlap.
  --run_ancestry            Specify whether to run ancestry.
  --add_sex                 Specify sample sex information
  
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

//Remove file extension, keeping path
def removeExtension(String path) {
    return path.replaceAll(/\.[^\.]+$/, '')
}

// Define workflow structure. Include some input/runtime tests here.
// See https://www.nextflow.io/docs/latest/dsl2.html?highlight=workflow#workflow
workflow {

// Show help message if --help is run or (||) a required parameter (input) is not provided

if ( params.help || !params.bamfile || !params.target_build || !params.ref || !params.dbsnp ||
    (params.target_build != 'GRCh37' && params.target_build != 'GRCh38') || 
    !(params.pgs_id || params.efo_id || params.pgp_id || params.scorefile) ) {  
// Invoke the help function above and exit
	switch (true) {
        case params.help:
            log.info "Help message:"

        case !params.bamfile:
            log.info "Error: Missing parameter '--bamfile'"

        case !params.target_build:
            log.info "Error: Missing parameter '--target_build'"

        case params.target_build != 'GRCh37' && params.target_build != 'GRCh38':
            log.info "Error: Invalid 'target_build'. Must be 'GRCh37' or 'GRCh38'"

        case !params.ref:
            log.info "Error: Missing parameter '--ref'"

        case !params.dbsnp:
            log.info "Error: Missing parameter '--dbsnp'"

        case !(params.pgs_id || params.efo_id || params.pgp_id || params.scorefile):
            log.info "Error: At least one of '--pgs_id', '--efo_id', '--pgp_id', or '--scorefile' must be provided"

        case !params.add_sex:
            log.info "WARNING: No sex data has been provided, if scorefile contains sex chromosomes pgsc_calc will fail loudly."

        default:
            // If none of the above are a problem, then run the workflow
            break
    }
    helpMessage()
    exit 1
	// consider adding some extra contigencies here.
	// could validate path of all input files in list?
	// could validate indexes for reference exist?

// If none of the above are a problem, then run the workflow
} else {
    if (!params.add_sex){
        log.info "WARNING: No sex data has been provided, if scorefile contains sex chromosomes pgsc_calc will fail loudly."
    }
	
	// DEFINE CHANNELS 
	// See https://www.nextflow.io/docs/latest/channel.html#channels
	// See https://training.nextflow.io/basic_training/channels/ 

	// DEMO CODE: DELETE FOR YOUR OWN WORKFLOWS - VALIDATE INPUT SAMPLES 
	//check_input(Channel.fromPath(params.input, checkIfExists: true))
	
	//SUBWORKFLOW Download - Scorefiles -----------------------------------------------
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

	//END SUBWORKFLOW -------------------------------------------------------------------

    //Reference File Check
    ref_file = file("${removeExtension(params.ref)}.fasta")
    if(!ref_file.exists()) {
        log.info "Reference file does not exist:"
        helpMessage()
        exit 1
    } else {
        ch_ref_fasta = Channel.fromPath(ref_file, checkIfExists: true)
    }

    //Reference index file check
    ref_idx_file = file("${removeExtension(params.ref)}.fasta.fai")
    if(!ref_idx_file.exists()) {
        // log.info "Reference fasta index file (.fasta.fai) does not exist, attempting to create it"
        //Create index file
        create_ref_idx("${removeExtension(params.ref)}.fasta")
        ch_ref_fasta_fai = create_ref_idx.out.ref_idx
    } else {
        ch_ref_fasta_fai = Channel.fromPath("${removeExtension(params.ref)}.fasta.fai", checkIfExists: true)
    }

    //Reference dictionary file check
    ref_dict_file = file("${removeExtension(params.ref)}.dict")
    if(!file(ref_dict_file).exists()) {
        // log.info "Reference dictionary file (.dict) does not exist, attempting to create it"
        //Create dict file
        create_ref_dict("${removeExtension(params.ref)}.fasta")
        ch_ref_dict = create_ref_dict.out.ref_dict
    } else {
        ch_ref_dict = Channel.fromPath(ref_dict_file, checkIfExists: true)
    }

    //Sample File Check
    bamfile_file = file(params.bamfile)
    if(!bamfile_file.exists()) {
        log.info "BAM file does not exist:"
        helpMessage()
        exit 1
    } else {
        ch_bamfile = Channel.fromPath(bamfile_file, checkIfExists: true)
    }

    //dbSNP File Check
    dbsnp_file = file(params.dbsnp)
    if (!dbsnp_file.exists()) {
        log.info "dbSNP file does not exist:"
        helpMessage()
        exit 1
    } else {
        ch_dbsnp = Channel.fromPath(dbsnp_file, checkIfExists: true)
    }
    
    //dbSNP Index File Check
    dbsnp_idx_file = file("${params.dbsnp}.idx")
    if (!dbsnp_idx_file.exists()) {
        // log.info "dbSNP index file does not exist, attempting to create it"
        //Create index file
        index_vcf(ch_dbsnp)
        ch_dbsnp_idx = index_vcf.out.vcf_idx
    } else {
        ch_dbsnp_idx = Channel.fromPath(dbsnp_idx_file, checkIfExists: true)
    }


	//Run generate_PRS_snp_positions_list for each scorefile
	generate_PRS_snp_positions_list(ch_scores.flatten(), ch_dbsnp)
		.collect()
		.set { prs_snp_positions_files }
	
    //Combine all PRS_snp_positions.list files into one, removing duplicate lines
	combine_PRS_snp_positions_lists(prs_snp_positions_files)

    //Remove 'chr' prefix from BAM file
    remove_chr_from_BAM(ch_bamfile)

	//Run GATK HaplotypeCaller on the BAM file
    GATK_haplotype_caller(
        remove_chr_from_BAM.out.bam_file_no_chr,
        ch_ref_fasta, 
        ch_ref_fasta_fai,
        ch_ref_dict,
        combine_PRS_snp_positions_lists.out.PRS_snp_positions,
        ch_dbsnp,
        ch_dbsnp_idx
        )

    //Run GATK GenotypeGVCFs on the gVCF file
    GATK_genotype_GVCFs(
        ch_ref_fasta,
        ch_ref_fasta_fai,
        ch_ref_dict,
        GATK_haplotype_caller.out.haplotypeCalled_gvcf,
        GATK_haplotype_caller.out.haplotypeCalled_gvcf_idx,
        ch_dbsnp,
        ch_dbsnp_idx,
        combine_PRS_snp_positions_lists.out.PRS_snp_positions
        )

    //Populate ALT Alleles, using custom logic (handing HOMREF, using dbSNP for indels)
    populate_alt_alleles(
        GATK_genotype_GVCFs.out.gvcf_genotyped,
        ch_scores.flatten().collect(),
        ch_ref_fasta,
        "${workflow.projectDir}/lib/populate_alt_alleles.sh",
        ch_dbsnp,
        combine_PRS_snp_positions_lists.out.PRS_snp_positions
    )

    //If sex information provided, convert VCF to PLINK format adding sex information
    if (params.add_sex) {
        //Normalise VCF to handle multiallelic sites
        bcftools_normalise_vcf(populate_alt_alleles.out.combined_processed_vcf, ch_ref_fasta)
        //Convert VCF to PLINK format, adding sex information
        convert_vcf_to_plink(bcftools_normalise_vcf.out.normalised_vcf, params.add_sex)
    } 
    //Create pgsc_calc samplesheet, either using PLINK or VCF
    make_samplesheet(
        params.add_sex ? convert_vcf_to_plink.out.pgen: populate_alt_alleles.out.combined_processed_vcf,
        params.add_sex,
        remove_chr_from_BAM.out.sample_name
    )

    //TODO -
    //Maybe I should run pgsc_calc test config with internet access, to load dependencies, then run main pipeline...

    // Create a channel of paths for the input files, either using PLINK or VCF files
    if(params.add_sex){
        populate_alt_alleles.out.combined_processed_vcf.concat(convert_vcf_to_plink.out.pgen, convert_vcf_to_plink.out.pvar, convert_vcf_to_plink.out.psam).set { ch_input_files }
    } else {
        populate_alt_alleles.out.combined_processed_vcf.set { ch_input_files }
    }
    
    //Run pg_sc_calc with input files
    pgsc_calc(
        make_samplesheet.out.samplesheet,
        ch_input_files.flatten().collect(),
        params.target_build,
        ch_scores.flatten().collect(),
        workflow.workDir,
        params.min_overlap,
        params.singularityCacheDir
    )
}}

// Print workflow execution summary 
workflow.onComplete {
    def traceFilePath = "${workflow.projectDir}/${params.outdir}/runInfo/trace-*.tsv"
    def command = "./lib/SU_utilisation_calculation.sh ${traceFilePath}"
    try {
        def process = command.execute()
        def su = process.text.trim()
        process.waitFor()
            summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
results     : ${params.outdir}
${su ? "\nSU utilisation: ${su}\n" : ""}
=======================================================================================
    """
    } catch (Exception e) {
        // If the command fails, print the summary without SU information
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
    }
    println summary
}
