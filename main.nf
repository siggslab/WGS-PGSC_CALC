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
include { check_input } from './modules/check_input'
include { group_samples } from './modules/group_samples'
include { generate_report } from './modules/generate_report' 

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
input       : ${params.input}
results     : ${params.outdir}
workDir     : ${workflow.workDir}
=======================================================================================

"""

/// Help function 
// This is an example of how to set out the help function that 
// will be run if run command is incorrect or missing. 

def helpMessage() {
    log.info"""
  Usage:  nextflow run main.nf --input <samples.tsv> 

  Required Arguments:

  --input		Specify full path and name of sample input file.

  Optional Arguments:

  --outdir	Specify path to output directory. 
	
""".stripIndent()
}

// Define workflow structure. Include some input/runtime tests here.
// See https://www.nextflow.io/docs/latest/dsl2.html?highlight=workflow#workflow
workflow {

// Show help message if --help is run or (||) a required parameter (input) is not provided

if ( params.help || params.input == false ){   
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
	check_input(Channel.fromPath(params.input, checkIfExists: true))

	// DEMO CODE: DELETE FOR YOUR OWN WORKFLOWS - EXAMPLE PROCESS - SPLIT SAMPLESHEET DEPENDING ON SEQUENCING PLATFORM
	// See https://training.nextflow.io/basic_training/processes/#inputs 
	// Define the input channel for this process
	group_samples_in = check_input.out.checked_samplesheet

	// Run the process with its input channel
	group_samples(group_samples_in)
	
	// DEMO CODE: DELETE FOR YOUR OWN WORKFLOWS - EXAMPLE PROCESS - SUMMARISE COHORT FROM SAMPLESHEETS
	// Define the input channel for this process using Nextflow mix operator and some groovy (the use of 'map')
	// See: https://www.nextflow.io/docs/latest/operator.html
	generate_report_in = group_samples.out.illumina
                     .map { file -> tuple(file, 'Illumina') }
                     .mix(group_samples.out.pacbio
                          .map { file -> tuple(file, 'PacBio') })
	
	// DEMO CODE: DELETE FOR YOUR OWN WORKFLOWS - Run the process with its input channel
	generate_report(generate_report_in)
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
