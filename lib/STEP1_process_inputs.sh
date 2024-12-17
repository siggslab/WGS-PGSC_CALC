#!/bin/bash

#rsync command to sync to nci
#rsync -avh --update /Users/bradford_el/Garvan_Internship/Michael_Procedure_handle_refonly eb9831@gadi.nci.org.au:/g/data/tn36/pgs/Michael_Procedure_handle_refonly; rsync -avh --update eb9831@gadi.nci.org.au:/g/data/tn36/pgs/Michael_Procedure_handle_refonly/Michael_Procedure_handle_refonly /Users/bradford_el/Garvan_Internship/


input_vcf=$1
scoring_file=$2
ref_file=$3

#Dependencies
module load bcftools
module load bedtools

#Seperate variant and reference sites from the input VCF file
bcftools view -v snps,indels,mnps,other -o variants.vcf $input_vcf
bcftools view -e 'TYPE="snp" | TYPE="indel" | TYPE="mnp" | TYPE="other"' -o non_variants.vcf $input_vcf

#Get input header
bcftools view -h $input_vcf > input_header.txt

# Extract the sample names from the VCF header and trim any trailing spaces/tabs
sample_name=$(bcftools query -l $input_vcf | xargs)

# Extract CHR, POS, REF, FORMAT, and SAMPLE fields from reference-only sites
bcftools query -f '%CHROM\t%POS\t%REF\t%FORMAT\t[%SAMPLE]\n' non_variants.vcf | \
awk -v sample="$sample_name" '
BEGIN {
    OFS="\t";
    # print "CHR", "POS", "REF", "FORMAT", sample
} 
{
    print $1, $2, $3, $4, $5
}' > reference_only_sites.tsv

#Extract harmonized CHR, POS, EFFECT, OTHER allele fields, and ID from the scoring file
zcat $scoring_file | awk '
BEGIN {
    FS = OFS = "\t";
}
# Skip lines starting with # or empty lines
/^#/ || /^$/ { next }
# Process the first non-comment line as the header
{
    if (header_found == 0) {
        for (i = 1; i <= NF; i++) {
            if ($i == "rsID") id_idx = i;
            else if ($i == "hm_chr") chr_idx = i;
            else if ($i == "hm_pos") pos_idx = i;
            else if ($i == "effect_allele") effect_allele_idx = i;
            else if ($i == "other_allele") other_allele_idx = i;
            else if ($i == "hm_inferOtherAllele" && !other_allele_idx) other_allele_idx = i;
        }
        # print "ID", "CHR", "POS", "EFFECT_ALLELE", "OTHER_ALLELE";
        header_found = 1;
        next;
    }
}
# Process the data lines
{
    id = (id_idx ? $id_idx : ".");
    print id, $chr_idx, $pos_idx, $effect_allele_idx, $other_allele_idx;
}' > scoring_file_extracted.tsv

# Sort both files by CHR and POS
sort -k1,1n -k2,2n reference_only_sites.tsv > reference_only_sites_sorted.tsv
sort -k2,2n -k3,3n scoring_file_extracted.tsv > scoring_file_extracted_sorted.tsv

# Merge the two tables by CHR and POS
awk '
BEGIN {
    FS = OFS = "\t";
    print "CHR", "POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE";
}
FNR==NR {
    key = $2 FS $3;
    data[key] = $0;
    next;
}
{
    key = $1 FS $2;
    if (key in data) {
        split(data[key], arr, FS);
        print $1, $2, arr[1], arr[4], arr[5], $3, $4, $5;
    }
}' scoring_file_extracted_sorted.tsv reference_only_sites_sorted.tsv > merged_table.tsv

# Split the merged table into SNPs and INDELs
awk '
BEGIN {
    FS = OFS = "\t";
    print "CHR", "POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE" > "merged_snps.tsv";
    print "CHR", "POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE" > "merged_indels.tsv";
}
NR > 1 {
    split($4, effect_alleles, "/");
    split($5, other_alleles, "/");
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
        print $0 > "merged_snps.tsv";
    } else {
        print $0 > "merged_indels.tsv";
    }
}' merged_table.tsv

# Function to get the reverse complement of a DNA sequence
# rev_comp() {
#     local seq=$1
#     echo "$seq" | tr 'ATCGatcg' 'TAGCtagc' | rev
# }

#TODO: Handle GT of './.'
# Assign ALT allele for SNPs
awk '
function rev_comp(seq, rev_seq) {
    rev_seq = "";
    for (i = length(seq); i > 0; i--) {
        base = substr(seq, i, 1);
        if (base == "A") rev_seq = rev_seq "T";
        else if (base == "T") rev_seq = rev_seq "A";
        else if (base == "C") rev_seq = rev_seq "G";
        else if (base == "G") rev_seq = rev_seq "C";
    }
    return rev_seq;
}
BEGIN {
    FS = OFS = "\t";
    print "CHR", "POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE", "ALT";
}
NR > 1 {
    ref = $6;
    effect = $4;
    other = $5;
    alt = "";

    if (effect == ref) {
        alt = other;
    } else if (other == ref) {
        alt = effect;
    } else {
        rev_comp_effect = rev_comp(effect);
        rev_comp_other = rev_comp(other);
        rev_comp_ref = rev_comp(ref);
        if (effect == rev_comp_ref) {
            alt = rev_comp_other;
        } else if (other == rev_comp_ref) {
            alt = rev_comp_effect;
        }
    }

    # Check for strand flip
    # if (alt!="" && (effect==alt && other==ref) ) {
    #     temp = ref;
    #     ref = alt;
    #     alt = temp;
    # }
    
    sample = $8;
    #Testing editing the GT field of the format
    # gt = substr(sample, 1, 3);
    # rest = substr(sample, 4);
    # sample = "1/1" rest;

    if (alt != "") {
        print $1, $2, $3, $4, $5, ref, $7, sample, alt;
        #print $1, $2, $3, $4, $5, ".", $7, $8, ".";
    } else {
        print "Discarding: Unable to determine ALT allele for SNP at " $1 ":" $2 > "/dev/stderr";
    }
}' merged_snps.tsv > merged_snps_with_alt.tsv

########SUBWORKFLOW: Add SEQ column to INDEL table #########
# Create a BED file from the INDELs table
awk '
BEGIN {
    FS = OFS = "\t";
}
NR > 1 {
    effect_len = length($4);
    other_len = length($5);
    max_len = (effect_len > other_len) ? effect_len : other_len;
    end_pos = $2 + max_len - 1;
    print $1, $2 - 1, end_pos, $1 ":" $2;
}' merged_indels.tsv > indels.bed

# Use bedtools to extract the genome sequences
bedtools getfasta -fi $ref_file -bed indels.bed -fo indels_sequences.fa

# Add the extracted sequences as a new column SEQ in the INDELs table
awk '
BEGIN {
    FS = OFS = "\t";
}
FNR == 1 && FILENAME == "merged_indels.tsv" {
    # Print the header line with the added SEQ column
    print $0, "SEQ";
    next;
}
NR == FNR {
    if ($0 ~ /^>/) {
        # Extract the key from the header line
        split(substr($0, 2), arr, "[:-]");
        key = arr[1] ":" arr[2] + 1;  # Adjust key to match CHR:POS format
    } else {
        # Store the sequence in the seq array
        seq[key] = $0;
    }
    next;
}
NR > 1 {
    # Construct the key using the chromosome and position
    key = $1 ":" $2;
    # Print the line with the added SEQ column
    print $0, (key in seq ? seq[key] : "");
}' indels_sequences.fa merged_indels.tsv > merged_indels_with_seq.tsv

#TODO: Handle REF with len > 1
#Assign ALT alleles for INDELs
awk '
BEGIN {
    FS = OFS = "\t";
    print "CHR", "POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE", "SEQ", "ALT";
}
NR > 1 {
    ref = $6;
    effect = $4;
    other = $5;
    alt = "";
    seq = $9;

    split(ref, ref_arr, "");
    split(effect, effect_arr, "");
    split(other, other_arr, "");

    if (effect_arr[1] != other_arr[1]) {
        print "Discarding: INDEL alleles do not match at " $1 ":" $2 ". Likely un-normalised INDEL" > "/dev/stderr";
        next;
    } else if ((effect_arr[1] != ref_arr[1]) && (other_arr[1] != ref_arr[1])) {
        print "Discarding: INDEL does not match on forward strand. Either on reverse strand or invalid at " $1 ":" $2 > "/dev/stderr";
        next;
    }

    seq_effect = substr(seq, 1, length(effect));
    seq_other = substr(seq, 1, length(other));

    if (seq_effect == effect && seq_other != other) {
        ref = effect;
        alt = other;
    } else if (seq_other == other && seq_effect != effect) {
        ref = other;
        alt = effect;
    } else if (seq_effect == effect && seq_other == other) {
        print "Discarding: Ambiguous INDEL at " $1 ":" $2 > "/dev/stderr";
        next;
    } else if (seq_effect != effect && seq_other != other) {
        print "Discarding: Unexpected invalid INDEL at " $1 ":" $2 > "/dev/stderr";
        next;
    }


    if (alt != "") {
        print $1, $2, $3, $4, $5, ref, $7, $8, $9, alt;
    } else {
        print "Discarding: Unable to determine ALT allele for INDEL at " $1 ":" $2 > "/dev/stderr";
    }
}' merged_indels_with_seq.tsv > merged_indels_with_seq_alt.tsv

##########SUBWORKFLOW - Compile ref only SNPs and INDELS into a VCF ######################
# Compile the VCF for SNPs
awk -v sample_name="$sample_name" '
BEGIN {
    FS = OFS = "\t";
    # print "##fileformat=VCFv4.2";
    # print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name;
    while ((getline line < "input_header.txt") > 0) {
        print line;
    }
}
NR > 1 {
    print $1, $2, $3, $6, $9, ".", "PASS", ".", $7, $8;
}' merged_snps_with_alt.tsv > reference_only_snps.vcf

# Compile the VCF for INDELs
awk -v sample_name="$sample_name" '
BEGIN {
    FS = OFS = "\t";
    # print "##fileformat=VCFv4.2";
    # print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name;
    while ((getline line < "input_header.txt") > 0) {
        print line;
    }
}
NR > 1 {
    print $1, $2, $3, $6, $10, ".", "PASS", ".", $7, $8;
}' merged_indels_with_seq_alt.tsv > reference_only_indels.vcf

#Compress and index VCF files
bcftools view -I reference_only_snps.vcf -O z -o reference_only_snps.vcf.gz
bcftools view -I reference_only_indels.vcf -O z -o reference_only_indels.vcf.gz
bcftools view -I variants.vcf -O z -o variants.vcf.gz
bcftools index reference_only_snps.vcf.gz
bcftools index reference_only_indels.vcf.gz
bcftools index variants.vcf.gz

# Combine the SNP and INDEL and variants VCF files
bcftools concat -a reference_only_snps.vcf.gz reference_only_indels.vcf.gz variants.vcf.gz -o combined_processed.vcf

