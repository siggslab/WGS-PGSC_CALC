#!/bin/bash

usage="bash benchmark.sh \\
    <sampleIDs> (a .txt file of newline seperated sampleIDs) \\
    <genotypingDir> (the path to the genotyping directory of bfiles (.bim, .bam, .fam)) \\
    <scorefile> (the path to the downloaded scorefile .txt.gz) \\
    <min_overlap> (a float 0-1 of the min_overlap to pass to pgsc_calc) \\
    <singularity_cache_dir> (the path to the cached singularity images) \\
    <wgs_data_dir> (the path to the wgs data, containing .bam files) \\
    <ref_path> (path to reference .fasta file) \\
    <sex_info_path> (path to .txt file with columns FID, IID, SEX) \\
    <dbsnp_path> (path to dbsnp file) \\
    <scorefile_arg> (scorefile_arg to pass to pgsc_calc eg. '--pgs_id PGS000001') \\
    <wgs_workflow_copy_dir> (path to local directory for WGS-PGSC_CALC Nextflow workflow) \\
    <genotype_workflow_copy_dir> (path to local directory for pgsc_calc Nextflow workflow)"

# Check if the required number of arguments is provided
if [ "$#" -ne 12 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $usage"
    echo "Missing parameters:"
    [ -z "$1" ] && echo "  <sampleIDs>"
    [ -z "$2" ] && echo "  <genotypingDir>"
    [ -z "$3" ] && echo "  <scorefile>"
    [ -z "$4" ] && echo "  <min_overlap>"
    [ -z "$5" ] && echo "  <singularity_cache_dir>"
    [ -z "$6" ] && echo "  <wgs_data_dir>"
    [ -z "$7" ] && echo "  <ref_path>"
    [ -z "$8" ] && echo "  <sex_info_path>"
    [ -z "$9" ] && echo "  <dbsnp_path>"
    [ -z "${10}" ] && echo "  <scorefile_arg>"
    [ -z "${11}" ] && echo "  <wgs_workflow_copy_dir>"
    [ -z "${12}" ] && echo "  <genotype_workflow_copy_dir>"
    exit 1
fi


#Dependencies
module load nextflow singularity
plink2_cmd="singularity exec singularity_images/plink2.sif plink2"

# Function to check if a file exists
check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 does not exist."
        exit 1
    fi
}

# Function to validate the samplesheet
validate_samplesheet() {
    samplesheet_path=$1

    # Check if the samplesheet path is correct
    check_file_exists $samplesheet_path

    # Read the header of the samplesheet
    header=$(head -n 1 $samplesheet_path)
    expected_header="sample,fq1,fq2,platform,library,center"

    # Check if the header matches the expected format
    if [ "$header" != "$expected_header" ]; then
        echo "Error: Invalid header format."
        exit 1
    fi

    # Check if all rows are CSV and not tab-delimited
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" == *$'\t'* ]]; then
            echo "Error: Found a column with a tab delimiter, expected CSV."
            exit 1
        fi
    done < <(tail -n +2 $samplesheet_path)  # Skip the header line
}

# Function to get sample IDs
get_sampleIDs() {
    sampleIDs_file=$1
    genotypingDir=$2
    output_file="sample_info.csv"

    while read line; do
        echo "Processing $line"
        for fam_file in ${genotypingDir}/PHASE_8.fam; do
            echo "Processing $fam_file"
            grep "$line" "$fam_file" | while read -r grep_line; do
                fam_file_base=$(basename "$fam_file" .fam)
                echo "$fam_file_base,$grep_line,$line" >> $output_file
            done
        done
    done < $sampleIDs_file
}

# Function to process PLINK files
plink_process() {
    sample_info_csv=$1
    genotypingDir=$2

    while IFS=, read -r bfile sample_fam sample_id; do
        echo "Processing bfile: $bfile, sample_fam: $sample_fam, sample_id: $sample_id"
        echo "$sample_fam" > temp_sample_fam.txt
        fam_id_from_fam=$(echo "$sample_fam" | awk '{print $1}')
        sample_id_from_fam=$(echo "$sample_fam" | awk '{print $2}')
        dir="${fam_id_from_fam}-${sample_id_from_fam}_pfiles"
        mkdir "$dir"
        $plink2_cmd -bfile "$genotypingDir/$bfile" --allow-extra-chr --chr 1-22,X,Y,XY -make-pgen --out "$dir/temp_${sample_id_from_fam}_pfile" --keep temp_sample_fam.txt 
        $plink2_cmd -pfile "$dir/temp_${sample_id_from_fam}_pfile" --allow-extra-chr --chr 1-22,X,Y,XY -make-pgen --out "$dir/${sample_id_from_fam}_pfile" --rm-dup force-first
        rm temp_sample_fam.txt
        rm $dir/temp*
    done < $sample_info_csv
}

# Function to calculate PGSC for genotyped samples
pgsc_calc_genotyped() {
    pfile_dirs=$1
    scorefile=$2
    min_overlap=$3
    singularity_cache_dir=$4
    main_nf=$5

    echo "sampleset,path_prefix,chrom,format" > samplesheet.csv
    for pfile_dir in ${pfile_dirs}; do
        pgen_file=$(ls ${pfile_dir} | grep ".pgen")
        pgen_basename=$(echo ${pgen_file} | sed 's/\.pgen$//')
        famid_iid=$(basename ${pfile_dir} | sed 's/_pfiles//')
        path_prefix="${pfile_dir}/${pgen_basename}"
        echo "$famid_iid,$path_prefix,,pfile" >> samplesheet.csv
    done

    export SINGULARITY_CACHEDIR=${singularity_cache_dir}
    export NXF_SINGULARITY_CACHEDIR=${singularity_cache_dir}
    export CONDA_PKGS_DIRS=${singularity_cache_dir}
    export HOME=$(pwd)/home
    mkdir -p $HOME

    nextflow run $main_nf \
        -profile gadi \
        --input samplesheet.csv \
        --target_build GRCh37 \
        --scorefile ${scorefile} \
        -w $(pwd)/work/pgsc_calc \
        --outDir $(pwd)/results \
        --min_overlap ${min_overlap} \
        --singularityCacheDir ${singularity_cache_dir}

    for file in results/*/score/aggregated_scores.txt.gz; do
        mv "$file" "${file%.txt.gz}_1.txt.gz"
    done
}
plink_process2() {
    sample_info_csv=$1
    imputed_pfile=$2

    while IFS=, read -r bfile sample_fam sample_id; do
        echo "Processing bfile: $bfile, sample_fam: $sample_fam, sample_id: $sample_id"
        echo "$sample_fam" > temp_sample_fam.txt
        fam_id_from_fam=$(echo "$sample_fam" | awk '{print $1}')
        sample_id_from_fam=$(echo "$sample_fam" | awk '{print $2}')
        dir="${fam_id_from_fam}-${sample_id_from_fam}_pfiles"
        mkdir "$dir"
        $plink2_cmd -bfile "$imputed_pfile" --allow-extra-chr --chr 1-22,X,Y,XY -make-pgen --out "$dir/${sample_id_from_fam}_pfile" --keep temp_sample_fam.txt 
        rm temp_sample_fam.txt
    done < $sample_info_csv
}

# Function to get WGS data paths
get_wgs_data_paths() {
    wgs_data_dir=$1
    output_file="bam_files.txt"

    find -L ${wgs_data_dir} -name "*.bam" > $output_file
    while read -r file; do
        ln -s "$file" .
    done < $output_file
}

# Function to run WGS PGSC calculation workflow
wgs_pgsc_calc_workflow() {
    ref_path=$1
    sex_info_path=$2
    singularity_cache_dir=$3
    dbsnp_path=$4
    scorefile_arg=$5
    wgs_workflow_main_nf=$6
    bam_path=$7
    workDir=$8
    min_overlap=$9

    export HOME=${workDir}/home
    export NXF_HOME=${workDir}/home
    mkdir -p $HOME
    rm -rf ${workDir}/work/wgs-pgsc-calc

    nextflow run ${wgs_workflow_main_nf} -profile gadi \
        ${scorefile_arg} \
        --bamfile "${bam_path}" \
        --target_build GRCh37 \
        --ref ${ref_path} \
        --singularityCacheDir ${singularity_cache_dir} \
        --dbsnp ${dbsnp_path} \
        --add_sex ${sex_info_path} \
        --min_overlap ${min_overlap} \
        -w ${workDir}/wgs-pgsc-calc \
        -resume
}

# Function to initiate results sheets
initiate_results_sheets() {
    echo "WGS RESULTS" > wgs_results.csv
    echo "FID,IID,PGS,SUM,DENOM,AVG,total_matched" >> wgs_results.csv
    echo "GENOTYPE RESULTS" > genotype_results.csv
    echo "FID,IID,PGS,SUM,DENOM,AVG,total_matched" >> genotype_results.csv
}

# Function to create results sheet
create_genotype_results_sheet() {
    genotype_match_file=$1
    genotype_score_file=$2
    genotype_results=$3

    zcat ${genotype_score_file} | tail -n +2 | while IFS=$'\t' read -r sampleset FID IID PGS SUM DENOM AVG; do
        total_matched=$(awk -F, -v fid="${FID}" -v iid="${IID}" '
            $4 == "matched" {
                sum += $12
            }
            END {
                print sum
            }' ${genotype_match_file})
        echo "${FID},${IID},${PGS},${SUM},${DENOM},${AVG},${total_matched}" >> ${genotype_results}
    done
}
create_wgs_results_sheet() {
    wgs_match_file=$1
    wgs_score_file=$2
    wgs_results=$3

    echo "Debug: wgs_match_file=${wgs_match_file}"
    echo "Debug: wgs_score_file=${wgs_score_file}"
    echo "Debug: wgs_results=${wgs_results}"

    zcat ${wgs_score_file} | tail -n +2 | while IFS=$'\t' read -r sampleset FID IID PGS SUM DENOM AVG; do
        total_matched=$(awk -F, -v fid="${FID}" -v iid="${IID}" '
            $4 == "matched" {
                sum += $12
            }
            END {
                print sum
            }' ${wgs_match_file})
        echo "Debug: FID=${FID}, IID=${IID}, PGS=${PGS}, SUM=${SUM}, DENOM=${DENOM}, AVG=${AVG}, total_matched=${total_matched}"
        echo "${FID},${IID},${PGS},${SUM},${DENOM},${AVG},${total_matched}" >> ${wgs_results}
    done
}

# Function to combine sheets
combine_sheets() {
    wgs_results=$1
    genotype_results=$2
    output_file="results.csv"

    cat $wgs_results > $output_file
    cat $genotype_results >> $output_file
}

# Main script execution ###############################################
sampleIDs=$1
genotypingDir=$2
scorefile=$3
min_overlap=$4
singularity_cache_dir=$5
wgs_data_dir=$6
ref_path=$7
sex_info_path=$8
dbsnp_path=$9
scorefile_arg=${10}
wgs_workflow_copy_dir=${11}
genotype_workflow_copy_dir=${12}


# Function to check if a file exists
check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 does not exist."
        exit 1
    fi
}

# Function to check if a directory exists
check_dir_exists() {
    if [ ! -d "$1" ]; then
        echo "Error: Directory $1 does not exist."
        exit 1
    fi
}

# Check if each required input is supplied and valid
check_file_exists "$sampleIDs"
check_dir_exists "$genotypingDir"
check_file_exists "$scorefile"
check_dir_exists "$singularity_cache_dir"
check_dir_exists "$wgs_data_dir"
check_file_exists "$ref_path"
check_file_exists "$sex_info_path"
check_file_exists "$dbsnp_path"
check_dir_exists "$wgs_workflow_copy_dir"
check_dir_exists "$genotype_workflow_copy_dir"



# Initiate results sheets
initiate_results_sheets

# Calculate PGSC for genotyped samples
j=1
for dir in ${genotypingDir}/*_pfiles; do
    temp_genotype_dir="temp_genotype_dir_${j}" && \
    mkdir -p $temp_genotype_dir/nextflow_files/ && \
    rsync -r --exclude 'singularity_cache/' $genotype_workflow_copy_dir $temp_genotype_dir/nextflow_files && \
    ln -sf ${genotype_workflow_copy_dir}/singularity_cache/ ${temp_genotype_dir}/nextflow_files/singularity_cache && \
    bash -c "
    pgsc_calc_genotyped() {
        pfile_dirs=\$1
        scorefile=\$2
        min_overlap=\$3
        singularity_cache_dir=\$4
        main_nf=\$5

        echo \"Debug: pfile_dirs=\${pfile_dirs}\"
        echo \"Debug: scorefile=\${scorefile}\"
        echo \"Debug: min_overlap=\${min_overlap}\"
        echo \"Debug: singularity_cache_dir=\${singularity_cache_dir}\"
        echo \"Debug: main_nf=\${main_nf}\"

        echo \"sampleset,path_prefix,chrom,format\" > samplesheet.csv
        
        for pfile_dir in \${pfile_dirs}; do
            echo \"Debug: Processing pfile_dir=\${pfile_dir}\">&2
            pgen_file=\$(ls \${pfile_dir} | grep \".pgen\")
            echo \"Debug: Found pgen_file=\${pgen_file}\">&2
            pgen_basename=\$(echo \${pgen_file} | sed 's/\.pgen$//')
            echo \"Debug: pgen_basename=\${pgen_basename}\">&2
            famid_iid=\$(basename \${pfile_dir} | sed 's/.*-\(.*\)_pfiles/\1/')
            echo \"Debug: famid_iid=\${famid_iid}\">&2
            path_prefix=\"\${pfile_dir}/\${pgen_basename}\"
            echo \"Debug: path_prefix=\${path_prefix}\">&2
            echo \"\$famid_iid,\$path_prefix,,pfile\" >> samplesheet.csv
        done

        export SINGULARITY_CACHEDIR=\${singularity_cache_dir}
        export NXF_SINGULARITY_CACHEDIR=\${singularity_cache_dir}
        export CONDA_PKGS_DIRS=\${singularity_cache_dir}
        export HOME=\$(pwd)/home
        mkdir -p \$HOME

        nextflow run \$main_nf \
            -profile gadi \
            --input samplesheet.csv \
            --target_build GRCh37 \
            --scorefile \${scorefile} \
            -w \$(pwd)/work/pgsc_calc \
            --outDir \$(pwd)/results \
            --min_overlap \${min_overlap} \
            --singularityCacheDir \${singularity_cache_dir} \
            --storage_account tn36,np30 \
            -resume

        #for file in results/*/score/aggregated_scores.txt.gz; do
        #    mv \"\$file\" \"\${file%.txt.gz}_1.txt.gz\"
        #done
        rm -rf results/results
        mkdir -p results/results
        mv results/\$famid_iid/* results/results/
    }
    cd $temp_genotype_dir && \
    echo \"Debug: dir=${dir} genotypingDir=${genotypingDir} scorefile=${scorefile} min_overlap=${min_overlap} temp_genotype_dir=${temp_genotype_dir} \" && \
    pgsc_calc_genotyped \"${dir}\" \"${scorefile}\" \"${min_overlap}\" \"/g/data/np30/users/eb9831/benchmark_wgs_np30/benchmarking_bash_script/${temp_genotype_dir}/nextflow_files/singularity_cache\" \"/g/data/np30/users/eb9831/benchmark_wgs_np30/benchmarking_bash_script/${temp_genotype_dir}/nextflow_files/main.nf\" && \
    cd .. " && \
    #Create ScoreFile
    create_genotype_results_sheet "${temp_genotype_dir}/results/results/match/*_summary.csv" "${temp_genotype_dir}/results/results/score/aggregated_scores.txt.gz" "genotype_results.csv" && \
    rm -rf $temp_genotype_dir/nextflow_files && \
    rm -rf $temp_genotype_dir/workflow_results && \
    mv $temp_genotype_dir/results/results/ $temp_genotype_dir/workflow_results && \
    rm -rf $temp_genotype_dir/results &
    j=$((j+1))
done

# Get WGS data paths
get_wgs_data_paths "$wgs_data_dir"

# Run WGS PGSC calculation workflow
i=1
while IFS= read -r bam_file; do
    temp_wgs_dir="temp_wgs_dir_${i}" && \
    rsync -r --exclude 'singularity_cache/' $wgs_workflow_copy_dir $temp_wgs_dir && \
    ln -sf ${wgs_workflow_copy_dir}/singularity_cache/ ${temp_wgs_dir}/WGS-PGSC_CALC/singularity_cache && \
    bash -c "
    wgs_pgsc_calc_workflow() {
        ref_path=\$1
        sex_info_path=\$2
        singularity_cache_dir=\$3
        dbsnp_path=\$4
        scorefile_arg=\$5
        wgs_workflow_main_nf=\$6
        bam_path=\$7
        workDir=\$8
        min_overlap=\$9

        export HOME=\${workDir}/home
        export NXF_HOME=\${workDir}/home
        mkdir -p \$HOME
        rm -rf \${workDir}/work/wgs-pgsc-calc

        nextflow run \${wgs_workflow_main_nf} -profile gadi \
            \${scorefile_arg} \
            --bamfile \"\${bam_path}\" \
            --target_build GRCh37 \
            --ref \${ref_path} \
            --singularityCacheDir \${singularity_cache_dir} \
            --dbsnp \${dbsnp_path} \
            --add_sex \${sex_info_path} \
            --min_overlap \${min_overlap} \
            -w \${workDir}/wgs-pgsc-calc \
            --storage_account tn36,np30 \
            -resume
    }

    cd $temp_wgs_dir && \
    wgs_pgsc_calc_workflow \"$ref_path\" \"$sex_info_path\" \"/g/data/np30/users/eb9831/benchmark_wgs_np30/benchmarking_bash_script/$temp_wgs_dir/WGS-PGSC_CALC/singularity_cache/\" \"$dbsnp_path\" \"$scorefile_arg\" \"WGS-PGSC_CALC/main.nf\" \"$bam_file\" \"$(pwd)\" \"$min_overlap\"&& \
    cd .." && \
    #Create ScoreFile
    create_wgs_results_sheet "${temp_wgs_dir}/results/results/*/match/*_summary.csv" "${temp_wgs_dir}/results/results/*/score/aggregated_scores.txt.gz" "wgs_results.csv" && \
    rm -rf $temp_wgs_dir/WGS-PGSC_CALC && \
    rm -rf $temp_wgs_dir/workflow_results && \
    mv $temp_wgs_dir/results/results/ $temp_wgs_dir/workflow_results && \
    rm -rf $temp_wgs_dir/results &
    i=$((i+1))
done < bam_files.txt

#Wait for all background jobs to finish
wait

# Combine sheets
combine_sheets "wgs_results.csv" "genotype_results.csv"

echo "Workflow completed successfully."