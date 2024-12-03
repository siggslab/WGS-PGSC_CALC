// Define the process
process combine_PRS_snp_positions_lists {
    // Define directives 
    // See: https://www.nextflow.io/docs/edge/process.html#processes
    debug = false //turn to true to print command stdout to screen
    publishDir "${params.outdir}/", mode: 'copy'
    container '' 

    // Define input 
    // See: https://www.nextflow.io/docs/latest/process.html#inputs
    input:
    path prs_snp_positions_files

    // Define output(s)
    // See: https://www.nextflow.io/docs/latest/process.html#outputs
    output:
    path("PRS_snp_positions.list"), emit: PRS_snp_positions

    // Define code to execute 
    // See: https://www.nextflow.io/docs/latest/process.html#script
    script:
    """
    # Combine all PRS_snp_positions files and remove duplicates
    cat ${prs_snp_positions_files} | sort | uniq > PRS_snp_positions.list

    # Debug: Print the content of the PRS_snp_positions.list file
    echo "Content of PRS_snp_positions.list:"
    cat PRS_snp_positions.list
    """
}