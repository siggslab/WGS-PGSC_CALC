# Pipeline name
WGS-PGSC_CALC
## Workflow description 
A pipeline for processing WGS data to calculate PRS using PGScatalog/pgsc_calc.

<img width="400" alt="Process Diagram" src="https://github.com/user-attachments/assets/0689c28e-3b62-46e5-9a99-2f1f3de2356f" />

## User guide 
Usage:  nextflow run main.nf --bamfile <bamfile.bam> --target_build <GRCh37|GRCh38> --ref <reference.fasta> --dbsnp <dbsnp.vcf> [--scorefile <scorefile.txt> | --pgs_id <pgs_id> | --efo_id <efo_id> | --pgp_id <pgp_id>] --outdir <output_directory>

    Required Arguments:
    
    --bamfile                 Specify full path and name of BAM file. 
    
    --target_build            Specify the target build. Must be 'GRCh37' or 'GRCh38'.
    
    --ref                     Specify full path and name of the reference genome file.
    
    --dbsnp                   Specify full path and name of the dbSNP file.
    
    --outdir                  Specify path to output directory.
    
    At least one of the following must be provided:
    
    --scorefile               Specify full path and name of score file.
    
    --pgs_id                  Specify PGS Catalog ID.
    
    --efo_id                  Specify EFO ID.
    
    --pgp_id                  Specify PGP ID.
     
    Optional Arguments:
    --help                    Show this help message and exit.
    
    --singularityCacheDir     Specify path to Singularity cache directory.
    
    --min_overlap             Specify minimum overlap.
    
    --run_ancestry            Specify whether to run ancestry.
    
    --add_sex                 Specify sample sex information
## Component tools 
PGScatalog/pgsc_calc

## Benchmark.sh Script
Benchmark.sh is located in lib/benchmark.sh, and is a script designed to benchmark this pipeline against pgsc_calc using imputed genotyping array data.
It requires locally cloning both this pipeline and the pgscatalog/pgsc_calc pipeline.

Script Usage: 

    bash benchmark.sh \

    <sampleIDs> (a .txt file of newline seperated sampleIDs) \
    
    <genotypingDir> (the path to the genotyping directory of bfiles (.bim, .bam, .fam)) \
    
    <scorefile> (the path to the downloaded scorefile .txt.gz) \
    
    <min_overlap> (a float 0-1 of the min_overlap to pass to pgsc_calc) \
    
    <singularity_cache_dir> (the path to the cached singularity images) \
    
    <wgs_data_dir> (the path to the wgs data, containing .bam files) \
    
    <ref_path> (path to reference .fasta file) \
    
    <sex_info_path> (path to .txt file with columns FID, IID, SEX) \
    
    <dbsnp_path> (path to dbsnp file) \
    
    <scorefile_arg> (scorefile_arg to pass to pgsc_calc eg. '--pgs_id PGS000001') \
    
    <wgs_workflow_copy_dir> (path to local directory for WGS-PGSC_CALC Nextflow workflow) \
    
    <genotype_workflow_copy_dir> (path to local directory for pgsc_calc Nextflow workflow)

## Acknowledgements/citations/credits
Lambert, S. A., Wingfield, B., Gibson, J. T., Gil, L., Ramachandran, S., Yvon, F., Saverimuttu, S., Tinsley, E., Lewis, E., Ritchie, S. C., Wu, J., Cánovas, R., McMahon, A., Harris, L. W., Parkinson, H., & Inouye, M. (2024). Enhancing the Polygenic Score Catalog with tools for score calculation and ancestry normalization. Nature Genetics, 56(10), 1989-1994. https://doi.org/10.1038/s41588-024-01937-x 
![image](https://github.com/user-attachments/assets/612dfb52-7a6f-40f5-9c5f-4606443540ee)
