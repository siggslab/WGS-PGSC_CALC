#!/bin/bash
# Bash script to process the scoring file and generate the PRS_snp_positions.list
scoring_file=$1
dbsnp_file=$2

if [[ -z "$scoring_file" || -z "$dbsnp_file" ]]; then
    echo "Usage: $0 <scoring_file> <dbsnp_file>"
    exit 1
fi

readCommand="cat"
if [[ $scoring_file == *.gz ]]; then
    readCommand="zcat"
fi

$readCommand ${scoring_file} | awk -v dbsnp_file="${dbsnp_file}" '
BEGIN {
    FS = "\t";
    OFS = "\n";
}
!/^#/ {
    if (header_found == 0) {
        header_found = 1;
        for (i=1; i<=NF; i++) {
            if ($i == "hm_rsID") hm_rsID_col = i;
            if ($i == "rsID") rsID_col = i;
        }
        next
    }
}
header_found == 1 {
    if (hm_rsID_col != "" && $hm_rsID_col != "") {
        rsIDs[NR] = $hm_rsID_col
    } else if (rsID_col != "" && $rsID_col != "") {
        rsIDs[NR] = $rsID_col
    }
}
END {
    for (id in rsIDs) {
        print rsIDs[id]
        #print "DEBUG: rsID in score file: " rsIDs[id] > "/dev/stderr"
    }
}' > temp_rsIDs


$readCommand ${scoring_file} | awk -v dbsnp_file="${dbsnp_file}" '
BEGIN {
    FS = "\t";
    # Read temp_rsIDs into an array
    while ((getline < "temp_rsIDs") > 0) {
        rsIDs[$0] = 1;
    }
    close("temp_rsIDs");
    
    while ((getline < dbsnp_file) > 0) {
        # Assuming dbsnp_file has fields: CHROM, POS, ID, REF, ALT
        if ($3 in rsIDs) {
            dbsnp[$3] = $1 ":" $2 "-" $2;  # ID -> CHROM:POS-POS mapping
        } 
    }
    close(dbsnp_file);
} 
!/^#/ { 
    if (header_found == 0) { 
        header_found = 1; 
        for (i=1; i<=NF; i++) { 
            if ($i == "hm_chr") chr_col = i; 
            if ($i == "hm_pos") pos_col = i; 
            if ($i == "hm_rsID") hm_rsID_col = i;
            if ($i == "rsID") rsID_col = i;
            if ($i == "other_allele") other_allele_col= i;
            else if ($i == "hm_inferOtherAllele" && !other_allele_col) other_allele_col = i;
            if ($i == "effect_allele") effect_allele_col = i;
        } 
        next 
    } 
} 
header_found == 1 { 
    split($effect_allele_col, effect_alleles, "/");
    split($other_allele_col, other_alleles, "/");
    is_snp = 1;
    for (i in effect_alleles) {
        if (length(effect_alleles[i]) > 1) {
            is_snp = 0;
            break;
        }
    }
    for (i in other_alleles) {
        if (length(other_alleles[i]) > 1) {
            is_snp = 0;
            break;
        }
    }
    if (is_snp) {
        if ($chr_col ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ && $pos_col ~ /^[0-9]+$/) {
            print $chr_col ":" $pos_col "-" $pos_col 
        }
    } else {
        # Determine the length of the indel
        #max_length = 0
        #for (i in effect_alleles) {
        #    if (length(effect_alleles[i]) > max_length) {
        #        max_length = length(effect_alleles[i])
        #    }
        #}
        #for (i in other_alleles) {
        #    if (length(other_alleles[i]) > max_length) {
        #        max_length = length(other_alleles[i])
        #    }
        #}
        #end_pos = $pos_col + max_length - 1
        #print "DEBUG: max_length=" max_length > "/dev/stderr"
        #print "DEBUG: end_pos=" end_pos > "/dev/stderr"
        #if ($chr_col ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ && $pos_col ~ /^[0-9]+$/) {
        #    print $chr_col ":" $pos_col "-" end_pos
        #}
        # Lookup chr and pos from dbsnp file using hm_rsID or rsID
        rsID = ($hm_rsID_col != "" ? $hm_rsID_col : $rsID_col)
        if (rsID != "") {
            if (rsID in dbsnp) {
                print "DEBUG: Using dbSNP for indel at: " dbsnp[rsID] ", rsID=" rsID > "/dev/stderr"
                print dbsnp[rsID]
            } else if ($chr_col ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ && $pos_col ~ /^[0-9]+$/) {
                print "rsID not in dbSNP file: " rsID > "/dev/stderr"
                print $chr_col ":" $pos_col "-" $pos_col 
            }
        } else if ($chr_col ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ && $pos_col ~ /^[0-9]+$/) {
            print "rsID not available for chr=" $chr_col ", pos=" $pos_col > "/dev/stderr"
            print $chr_col ":" $pos_col "-" $pos_col 
        }

    }
}' > $(basename "${scoring_file}" .txt.gz | sed 's/\.txt$//')_PRS_snp_positions.list
