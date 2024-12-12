// Define the process
process correct_alt_for_homref {
	// Define directives 
	// See: https://nextflow.io/docs/edge/process.html#processes
	debug = false //turn to true to print command stdout to screen
	tag "" 
	publishDir "${params.outdir}/", mode: 'copy'
  container '' 

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path(input_vcf)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path("${input_vcf.baseName}_homref_corrected.vcf"), emit: homref_corrected_vcf

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	awk '
        BEGIN {
            OFS = "\\t"  # Set the output field separator to a tab
        }
        /^#/ {
            # Print header lines as-is
            print
            next
        }
        {
            # Split the FORMAT column into fields
            split(\$9, format_fields, ":")
            
            # Find the index of the GT field
            gt_index = -1
            for (i = 1; i <= length(format_fields); i++) {
                if (format_fields[i] == "GT") {
                    gt_index = i
                    break
                }
            }
            
            # If GT field is found, process sample columns
            if (gt_index != -1) {
                # Split the first sample column (column 10) into fields
                split(\$10, sample_fields, ":")
                
                # Check if GT is 0/0 and REF is not "."
                if (sample_fields[gt_index] == "0/0" && \$4 != ".") {
                    #Set ALT to the other three bases
					#Solution proposed here: https://github.com/PGScatalog/pgsc_calc/discussions/247#discussioncomment-8596664
					#\$5 = \$4
					if (\$4 == "A") {
						\$5 = "T,C,G"
					} else if (\$4 == "T") {
						\$5 = "A,C,G"
					} else if (\$4 == "C") {
						\$5 = "A,T,G"
					} else if (\$4 == "G") {
						\$5 = "A,T,C"
					}
					#if (\$4 == "A") {
					#	\$5 = "C"
					#} else if (\$4 == "T") {
					#	\$5 = "A"
					#} else if (\$4 == "C") {
					#	\$5 = "G"
					#} else if (\$4 == "G") {
					#	\$5 = "C"
					#}
                }
            }
            
            # Print the modified or unmodified line
            print
        }
    ' "${input_vcf}" > "${input_vcf.baseName}_homref_corrected.vcf"
	"""
 }