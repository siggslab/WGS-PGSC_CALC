#!/bin/bash

input_vcf=$1
scoring_file=$2
ref_file=$3
dbsnp_file=$4
snp_list_file=$5

#Initialise Logging Variables
num_total=0
num_variants=0
num_non_variants=0
num_non_variants_not_in_scoring_file=0
num_non_variant_snps=0
num_non_variant_extra_multiallelic=0
num_non_variant_indels=0
num_non_variant_snps_assigned_alts=0
num_non_variant_indels_assigned_alts=0
num_non_variant_discarded_snps=0
num_non_variant_discarded_indels=0

#Dependencies
module load bcftools
module load bedtools

#Count total number of sites in input VCF
num_total=$(bcftools view -H $input_vcf | wc -l)

#Seperate variant and reference sites from the input VCF file
bcftools view -v snps,indels,mnps,other -o variants.vcf $input_vcf
#-T ^variants.vcf ensures no sites in variants.vcf are included in non_variants.vcf ( some multiallelics were origninally being included)
bcftools view -e 'TYPE="snp" | TYPE="indel" | TYPE="mnp" | TYPE="other"' -T ^variants.vcf -o non_variants.vcf $input_vcf 

#Count number of variants and non-variants
num_variants=$(bcftools view -H variants.vcf | sort | uniq | wc -l)
num_non_variants=$(bcftools view -H non_variants.vcf | sort | uniq |  wc -l)

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
echo "Scoring files: $scoring_file"
for file in $scoring_file; do
    echo "Processing $file"
    if [[ $file == *.gz ]]; then
        readcommand="zcat"
    elif [[ $file == *.txt ]]; then
        readcommand="cat"
    fi
    $readcommand $file | awk '
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
                else if ($i == "hm_chr") hm_chr_idx = i;
                else if ($i == "hm_pos") hm_pos_idx = i;
                else if ($i == "chr_name") chr_idx = i;
                else if ($i == "chr_position") pos_idx = i;
                else if ($i == "effect_allele") effect_allele_idx = i;
                else if ($i == "other_allele") other_allele_idx = i;
                else if ($i == "hm_inferOtherAllele" && !other_allele_idx) other_allele_idx = i;
            }
            # print "ID", "HM_CHR", "HM_POS", "CHR", "POS", "EFFECT_ALLELE", "OTHER_ALLELE";
            header_found = 1;
            print "id_idx=" id_idx, "hm_chr_idx=" hm_chr_idx, "hm_pos_idx=" hm_pos_idx, "chr_idx=" chr_idx, "pos_idx=" pos_idx, "effect_allele_idx=" effect_allele_idx, "other_allele_idx=" other_allele_idx > "/dev/stderr";
            next;
        }
    }
    # Process the data lines
    {
        if(header_found){
            id = (id_idx ? $id_idx : ".");
            hm_chr = (hm_chr_idx ? $hm_chr_idx : ".");
            hm_pos = (hm_pos_idx ? $hm_pos_idx : ".");
            chr = (chr_idx ? $chr_idx : ".");
            pos = (pos_idx ? $pos_idx : ".");

            print id, hm_chr, hm_pos, chr, pos, $effect_allele_idx, $other_allele_idx;
        }
    }' >> scoring_file_extracted.tsv
done

# Remove duplicate lines from scoring_file_extracted.tsv
awk '!seen[$0]++' scoring_file_extracted.tsv > scoring_file_extracted_unique.tsv
mv scoring_file_extracted_unique.tsv scoring_file_extracted.tsv


# Sort both files by CHR and POS
sort -k1,1n -k2,2n reference_only_sites.tsv > reference_only_sites_sorted.tsv
sort -k2,2n -k3,3n scoring_file_extracted.tsv > scoring_file_extracted_sorted.tsv

# Merge the two tables by CHR and POS
awk -v num_non_variants_not_in_scoring_file="$num_non_variants_not_in_scoring_file" '
BEGIN {
    FS = OFS = "\t";
    print "HM_CHR", "HM_POS", "score_CHR", "score_POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE";
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
        #print "Debug: $1=" $1 ", $2=" $2 ", $3=" $3 ", $4=" $4 ", $5=" $5 ", $6=" $6 ", $7=" $7 ", arr[1]=" arr[1] ", arr[2]=" arr[2] ", arr[3]=" arr[3] ", arr[4]=" arr[4] ", arr[5]=" arr[5] > "/dev/stderr";
        #ARR is scoring_file
        #$ is reference_only_sites
        #ref_only format: "CHR", "POS", "REF", "FORMAT", sample
        #score format: "ID", "HM_CHR", "HM_POS", "CHR", "POS", "EFFECT_ALLELE", "OTHER_ALLELE"
        
        print $1, $2, arr[4], arr[5], arr[1], arr[6], arr[7], $3, $4, $5;
    } else {
        num_non_variants_not_in_scoring_file++;
    }
}
END{
    print num_non_variants_not_in_scoring_file > "num_non_variants_not_in_scoring_file.tmp";
}' scoring_file_extracted_sorted.tsv reference_only_sites_sorted.tsv > merged_table.tsv

# Read the number of non_variants not in the scoring file
num_non_variants_not_in_scoring_file=$(cat num_non_variants_not_in_scoring_file.tmp)
rm num_non_variants_not_in_scoring_file.tmp

# Split the merged table into SNPs and INDELs
awk -v num_non_variant_snps="$num_non_variant_snps" -v num_non_variant_indels="$num_non_variant_indels" -v num_non_variant_extra_multiallelic="$num_non_variant_extra_multiallelic" '
BEGIN {
    FS = OFS = "\t";
    print "HM_CHR", "HM_POS", "score_CHR", "score_POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE" > "merged_snps.tsv";
    print "HM_CHR", "HM_POS", "score_CHR", "score_POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE" > "merged_indels.tsv";
}
NR > 1 {
    split($6, effect_alleles, "/");
    split($7, other_alleles, "/");
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
    # for (i in other_alleles){
    #     if (i != 1) {
    #         num_non_variant_extra_multiallelic++;
    #     }
        
        if (is_snp) {
            # print $1, $2, $3, $4, $5, effect_alleles[1], other_alleles[i], $8, $9, $10 > "merged_snps.tsv";
            num_non_variant_snps++;
            print $0 > "merged_snps.tsv";
        } else {
            # print $1, $2, $3, $4, $5, effect_alleles[1], other_alleles[i], $8, $9, $10 > "merged_indels.tsv";
            num_non_variant_indels++;
            print $0 > "merged_indels.tsv";
        }
    # }
}
END {
    print num_non_variant_snps > "num_non_variant_snps.tmp";
    print num_non_variant_indels > "num_non_variant_indels.tmp";
    print num_non_variant_extra_multiallelic > "num_non_variant_extra_multiallelic.tmp";
}' merged_table.tsv

# Read the updated values from the temporary files
num_non_variant_snps=$(cat num_non_variant_snps.tmp)
num_non_variant_indels=$(cat num_non_variant_indels.tmp)
num_non_variant_extra_multiallelic=$(cat num_non_variant_extra_multiallelic.tmp)

# Remove the temporary files
rm num_non_variant_snps.tmp num_non_variant_indels.tmp num_non_variant_extra_multiallelic.tmp

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
    print "HM_CHR", "HM_POS", "score_CHR", "score_POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE", "ALT";
    num_non_variant_snps_assigned_alts = 0;
    num_non_variant_discarded_snps = 0;
}
NR > 1 {
    ref = $8;
    effect = $6;
    other = $7;
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

    sample = $10;

    if (alt != "") {
        print $1, $2, $3, $4, $5, $6, $7, ref, $9, sample, alt;
        num_non_variant_snps_assigned_alts++;
    } else {
        print "Discarding: Unable to determine ALT allele for SNP at " $1 ":" $2 > "/dev/stderr";
        num_non_variant_discarded_snps++;
    }
}
END {
    print num_non_variant_snps_assigned_alts > "num_non_variant_snps_assigned_alts.tmp";
    print num_non_variant_discarded_snps > "num_non_variant_discarded_snps.tmp";
}' merged_snps.tsv > merged_snps_with_alt.tsv

# Read the updated values from the temporary files
num_non_variant_snps_assigned_alts=$(cat num_non_variant_snps_assigned_alts.tmp)
num_non_variant_discarded_snps=$(cat num_non_variant_discarded_snps.tmp)

# Remove the temporary files
rm num_non_variant_snps_assigned_alts.tmp num_non_variant_discarded_snps.tmp

########SUBWORKFLOW: Add SEQ column to INDEL table #########
# Create a BED file from the INDELs table
awk '
BEGIN {
    FS = OFS = "\t";
}
NR > 1 {
    effect_len = length($6);
    other_len = length($7);
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
        #TODO: check if the +1 is to do with harmonised/not
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


# Check if at least one row has a value for id
if awk 'NR > 1 && $5 != "." { found=1; exit } END { exit !found }' merged_indels_with_seq.tsv; then
    echo "At least one row of input indels has an rsID"
    # First, create a temporary file with the header
    awk '
    BEGIN {
        FS = OFS = "\t";
    }
    FNR == 1 && FILENAME == ARGV[2] {
        # Print the header line with the added dbsnp_alt column
        print $0, "dbsnp_alt";
        next;
    }' $dbsnp_file merged_indels_with_seq.tsv > merged_indels_with_seq_dbsnp_alt.tsv
    
    # Load SNP list file into an array
    declare -A snp_list
    while IFS= read -r line; do
        snp_list["$line"]=1
    done < "$snp_list_file"
    
    # Process the dbsnp file and store the ALT values in a temporary file
    awk -v snp_list_file="$snp_list_file" '
    BEGIN {
        FS = OFS = "\t";
        while ((getline line < snp_list_file) > 0) {
            snp_list[line] = 1;
        }
        close(snp_list_file);
    }
    FNR > 1 && $0 !~ /^#/ {
        # Data is in the format: CHROM, POS, REF, ALT
        data = $1 ":" $2 ":" $4 ":" $5;
        # Construct the snp_list_key using the chromosome and position range
        snp_list_key = $1 ":" $2 "-" $2;
        # Print the key and ALT value if snp_list_key is in snp_list
        if (snp_list_key in snp_list) {
            print $3, data;
        }
    }' $dbsnp_file > dbsnp_temp.tsv

    # Process the merged_indels_with_seq.tsv file and add the dbsnp_alt column based on CHR, POS columns
    awk '
    BEGIN {
        FS = OFS = "\t";
        print "Debug: Starting processing indels" > "/dev/stderr";
    }
    FILENAME == "dbsnp_temp.tsv" {
        # Read the dbsnp_temp.tsv file and store the ALT values in an array
        dbsnp_data[$1] = $2;
        next;
    }
    FILENAME == "merged_indels_with_seq.tsv" && FNR > 1 {
        #FORMAT HM_CHR  HM_POS  score_CHR  score_POS rsID EFFECT_ALLELE OTHER_ALLELE REF FORMAT SAMPLE SEQ
        
        # Get CHROM, POS, REF, ALT from rsID
        data = dbsnp_data[$5];

        if (data != "") {
            split(data, arr, ":");
            dbsnp_chr = arr[1];
            dbsnp_pos = arr[2];
            dbsnp_ref = arr[3];
            dbsnp_alt = arr[4];
            print "Debug: dbsnp_chr=" dbsnp_chr ", dbsnp_pos=" dbsnp_pos ", dbsnp_ref=" dbsnp_ref ", dbsnp_alt=" dbsnp_alt > "/dev/stderr";
            if( ((dbsnp_chr == $1) && (dbsnp_pos == $2) ) || ( (dbsnp_chr == $3) && (dbsnp_pos == $4)) ){
                print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, dbsnp_alt;
            } else {
                print "Debug: dbsnp pos does not match scorefile." > "/dev/stderr";
                print $0, ".";
            }
        } else {
            print "Debug: rsID " $5 " not found in dbsnp file" > "/dev/stderr";
            print $0, ".";
        }
        
    }' dbsnp_temp.tsv merged_indels_with_seq.tsv >> merged_indels_with_seq_dbsnp_alt.tsv

    # Clean up temporary files
    rm dbsnp_temp.tsv
else
    echo "No rows of input indels have an rsID"
    
    #Add dbSNP header to file
    awk '
    BEGIN {
        FS = OFS = "\t";
    }
    FNR == 1  {
        # Print the header line with the added dbsnp_alt column
        print $0, "dbsnp_alt";
        next;
    }
    FNR > 1{print $0, "."}
    ' merged_indels_with_seq.tsv > merged_indels_with_seq_dbsnp_alt.tsv
fi

#####Add alt based on dbsnp_alt column
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
    print "HM_CHR", "HM_POS", "score_CHR", "score_POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE", "SEQ", "dbsnp_alt", "ALT";
    num_non_variant_indels_assigned_alts = 0;
    num_non_variant_discarded_indels = 0;
}
NR > 1 && $12 != "." {
    ref = $8;
    effect = $6;
    other = $7;
    alt = "";
    seq = $11;
    dbsnp_alt = $12;

    split(ref, ref_arr, "");
    split(effect, effect_arr, "");
    split(other, other_arr, "");

    split(dbsnp_alt, dbsnp_alt_arr, ",");

    #Peform sanity check to ensure dbsnp alt is ok to use
    ok_flag = 0;
    if ((ref == other && dbsnp_alt == effect) || 
        (ref == effect && dbsnp_alt == other) || 
        (ref == rev_comp(other) && dbsnp_alt == rev_comp(effect)) || 
        (ref == rev_comp(effect) && adbsnp_alt == rev_comp(other))) {
        ok_flag = 1;
        alt=dbsnp_alt;
    } else {
        for (i in dbsnp_alt_arr) {
            if ((ref == other && dbsnp_alt_arr[i] == effect) || 
            (ref == effect && dbsnp_alt_arr[i] == other) || 
            (ref == rev_comp(other) && dbsnp_alt_arr[i] == rev_comp(effect)) || 
            (ref == rev_comp(effect) && dbsnp_alt_arr[i] == rev_comp(other))) {
                ok_flag = 1;
                alt = dbsnp_alt_arr[i];
                break;
            }
        }
    }

    if(ok_flag){
        print $1, $2, $3, $4, $5, $6, $7, ref, $9, $10, $11, dbsnp_alt, alt;
        num_non_variant_indels_assigned_alts++;
    } else {
        dbsnp_alt = ".";
        print $1, $2, $3, $4, $5, $6, $7, ref, $9, $10, $11, dbsnp_alt, alt;
    }
    
}
NR > 1 && $12 == "."{print $0;}
END {
    print num_non_variant_indels_assigned_alts > "num_non_variant_indels_assigned_alts.tmp";
    print num_non_variant_discarded_indels > "num_non_variant_discarded_indels.tmp";
}' merged_indels_with_seq_dbsnp_alt.tsv > merged_indels_with_seq_dbsnp_alt_processed.tsv

# Read the updated logging values from the temporary files
num_non_variant_indels_assigned_alts=$((num_non_variant_indels_assigned_alts + $(cat num_non_variant_indels_assigned_alts.tmp)))
num_non_variant_discarded_indels=$((num_non_variant_discarded_indels + $(cat num_non_variant_discarded_indels.tmp)))

# Remove the temporary files
rm num_non_variant_indels_assigned_alts.tmp num_non_variant_discarded_indels.tmp

#TODO: Handle REF with len > 1
# Assign ALT alleles for INDELs, where ALT is not pulled from dbSNP
awk '
BEGIN {
    FS = OFS = "\t";
    print "HM_CHR", "HM_POS", "score_CHR", "score_POS", "rsID", "EFFECT_ALLELE", "OTHER_ALLELE", "REF", "FORMAT", "SAMPLE", "SEQ", "dbsnp_alt", "ALT";
    num_non_variant_indels_assigned_alts = 0;
    num_non_variant_discarded_indels = 0;
}
NR > 1 && $12 == "." {
    ref = $8;
    effect = $6;
    other = $7;
    alt = "";
    seq = $11;

    split(ref, ref_arr, "");
    split(effect, effect_arr, "");
    split(other, other_arr, "");

    if (effect_arr[1] != other_arr[1]) {
        print "Discarding: INDEL alleles do not match at " $1 ":" $2 ". Likely un-normalised INDEL" > "/dev/stderr";
        num_non_variant_discarded_indels++;
        next;
    } else if ((effect_arr[1] != ref_arr[1]) && (other_arr[1] != ref_arr[1])) {
        print "Discarding: INDEL does not match on forward strand. Either on reverse strand or invalid at " $1 ":" $2 > "/dev/stderr";
        num_non_variant_discarded_indels++;
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
        num_non_variant_discarded_indels++;
        next;
    } else if (seq_effect != effect && seq_other != other) {
        print "Discarding: Unexpected invalid INDEL at " $1 ":" $2 > "/dev/stderr";
        num_non_variant_discarded_indels++;
        next;
    }

    if (alt != "") {
        print $1, $2, $3, $4, $5, $6, $7, ref, $9, $10, $11, $12, alt;
        num_non_variant_indels_assigned_alts++;
    } else {
        print "Discarding: Unable to determine ALT allele for INDEL at " $1 ":" $2 > "/dev/stderr";
        num_non_variant_discarded_indels++;
    }
}
NR > 1 && $12 != "."{
    print $0;
}
END {
    print num_non_variant_indels_assigned_alts > "num_non_variant_indels_assigned_alts.tmp";
    print num_non_variant_discarded_indels > "num_non_variant_discarded_indels.tmp";
}' merged_indels_with_seq_dbsnp_alt_processed.tsv > merged_indels_with_seq_dbsnp_alt_processed_alt.tsv

# Read the updated values from the temporary files
num_non_variant_indels_assigned_alts=$((num_non_variant_indels_assigned_alts + $(cat num_non_variant_indels_assigned_alts.tmp)))
num_non_variant_discarded_indels=$((num_non_variant_discarded_indels + $(cat num_non_variant_discarded_indels.tmp)))

# Remove the temporary files
rm num_non_variant_indels_assigned_alts.tmp num_non_variant_discarded_indels.tmp

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
    print $1, $2, $5, $8, $11, ".", "PASS", ".", $9, $10;
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
    print $1, $2, $5, $8, $13, ".", "PASS", ".", $9, $10;
}' merged_indels_with_seq_dbsnp_alt_processed_alt.tsv > reference_only_indels.vcf

#Compress and index VCF files
bcftools view -I reference_only_snps.vcf -O z -o reference_only_snps.vcf.gz
bcftools view -I reference_only_indels.vcf -O z -o reference_only_indels.vcf.gz
bcftools view -I variants.vcf -O z -o variants.vcf.gz
bcftools index reference_only_snps.vcf.gz
bcftools index reference_only_indels.vcf.gz
bcftools index variants.vcf.gz

# Combine the SNP and INDEL and variants VCF files
bcftools concat -a reference_only_snps.vcf.gz reference_only_indels.vcf.gz variants.vcf.gz -o combined_processed.vcf


# Print the logs
{
echo "Total number of variants in input VCF: $num_total"
echo "Number of variant sites: $num_variants"
echo "Number of non-variant sites: $num_non_variants"
echo "Number of non-variant sites in scoring file: $(($num_non_variants - $num_non_variants_not_in_scoring_file))"
echo "Number of non-variant SNPs: $num_non_variant_snps"
if [ $num_non_variant_extra_multiallelic -gt 0 ]; then
    echo "$(($num_non_variant_snps - $num_non_variant_extra_multiallelic)) plus $num_non_variant_extra_multiallelic expanded multiallelics"
fi
echo "Number of non-variant INDELs: $num_non_variant_indels"
echo "Number of non-variant SNPs assigned ALT alleles: $num_non_variant_snps_assigned_alts"
echo "Number of non-variant SNPs discarded: $num_non_variant_discarded_snps"
echo "Number of non-variant INDELs assigned ALT alleles: $num_non_variant_indels_assigned_alts"
echo "Number of non-variant INDELs discarded: $num_non_variant_discarded_indels"
} | tee populate_alt_alleles.log

