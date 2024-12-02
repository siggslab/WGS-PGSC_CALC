process check_input {
    tag "INPUT: ${input.fileName}"

    input:
		path(input)
    
    output:
    path("validated_samplesheet.csv") , emit: checked_samplesheet

    script: // This process runs ../bin/samplesheetchecker.sh
    // See example at: https://github.com/Sydney-Informatics-Hub/Parabricks-Genomics-nf/blob/main/bin/samplesheetchecker.py
    """
    samplesheetchecker.sh \\
      ${input} \\
      > validated_samplesheet.txt
    """
}