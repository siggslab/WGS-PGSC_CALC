// Define the process
process generate_PRS_snp_positions_list {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "INPUT: ${scoring_file}" 
	publishDir "${params.outdir}/", mode: 'copy'
    container '' 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(scoring_file)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("${scoring_file.simpleName}_PRS_snp_positions.list"), emit: PRS_snp_positions

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	// Determine the command to read the file based on its extension
	def readCommand = scoring_file.name.endsWith(".gz") ? "zcat" : "cat"
	"""
	# Bash script to process the scoring file and generate the PRS_snp_positions.list
	$readCommand ${scoring_file} | awk 'BEGIN {FS="\t"} 
		!/^#/ { 
			if (header_found == 0) { 
				header_found = 1; 
				for (i=1; i<=NF; i++) { 
					if (\$i == "hm_chr") chr_col = i; 
					if (\$i == "hm_pos") pos_col = i; 
				} 
				next 
			} 
		} 
		header_found == 1 { 
			print \$chr_col \":\" \$pos_col \"-\" \$pos_col 
		}' > ${scoring_file.simpleName}_PRS_snp_positions.list
	"""
 }