process generate_report {
    tag "INPUT: ${cohort}"
		publishDir "${params.outdir}", mode: 'copy'

	  input:
    tuple path(samplesheet), val(cohort)

    output:
    path("${cohort}_summary.txt")

    script:
    """
    # Count the number of samples in the cohort
    total_samples=\$(awk 'NF>0 {count++} END {print count+0}' ${samplesheet})

    # Count the number of fastq pairs for each sample
    total_fq_pairs=\$(awk -F, 'NF>0 {if(\$2 != "" && \$3 != "") count++} END {print count+0}' ${samplesheet})

    # Generate the cohort summary report
    echo "${cohort} Cohort Summary" > ${cohort}_summary.txt
    echo "===========================" >> ${cohort}_summary.txt
    echo "Total samples: \$total_samples" >> ${cohort}_summary.txt
    echo "Total fastq pairs per sample: \$total_fq_pairs" >> ${cohort}_summary.txt

    # List the files associated with each sample
    echo "" >> ${cohort}_summary.txt
    echo "Sample details:" >> ${cohort}_summary.txt
    awk -F, '{print "Sample: "\$1", Fastq 1: "\$2", Fastq 2: "\$3}' ${samplesheet} >> ${cohort}_summary.txt
    """
}