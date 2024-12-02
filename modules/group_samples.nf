process group_samples {
    tag "INPUT: ${checked_samplesheet.fileName}"
	publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(checked_samplesheet)

    output:
    path("samplesheet_illumina.csv"), emit: illumina
    path("samplesheet_pacbio.csv"), emit: pacbio

    script:
    """
    awk -F, '\$4 == "illumina"' OFS=, ${checked_samplesheet} > samplesheet_illumina.csv
    awk -F, '\$4 == "pacbio"' OFS=, ${checked_samplesheet} > samplesheet_pacbio.csv
    """
}
