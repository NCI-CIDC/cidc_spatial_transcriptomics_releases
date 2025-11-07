#!/usr/bin/env python3

## Generate the `--sample` parameter for `spaceranger count` based on the number of unique sample names
## present in the FASTQ directory. Multiple sample names indicate sequencing across multiple flow cells.
## FASTQ file naming convention for Space Ranger input
## Refer to https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/fastqs-specifying-fastqs#nonstandard-names
## [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz where Read Type is one of:
## I1: Sample index read (optional)
## I2: Sample index read (optional)
## R1: Read 1
## R2: Read 2


## Import modules
import glob
import os
import re

## Log output and stderr
log_file = open(snakemake.log[0], "w")
sys.stdout = log_file
sys.stderr = log_file

## Directory location that contains the FASTQs
fastqs_dir = snakemake.input.fastqs_dir
print(fastqs_dir)

## Use glob to get all .fastq.gz files in the directory
fastq_files = glob.glob(os.path.join(fastqs_dir, '*.fastq.gz'))

## Extract the base file name for the FASTQs
fastq_names = [os.path.basename(file) for file in fastq_files]

## Print the list of files
for fastq_name in fastq_names:
    print(fastq_name)

## Extract sample names dynamically
sample_names = set()
for fastq_name in fastq_names:
    sample_name = re.sub(r'_S\d+_.*', '', fastq_name)  ## Removes _S<number>_ dynamically
    sample_names.add(sample_name)

unique_sample_names = (sorted(sample_names))

## If multiple sample names exist, join them into a comma-separated string.
sample_name_str = "".join(unique_sample_names) if len(unique_sample_names) == 1 else ",".join(unique_sample_names)
print(sample_name_str)

if len(unique_sample_names) > 1:
    print("There is more than one sample name for the FASTQs.")
else:
    print("There is only one sample name for the FASTQs.")

## Write the sample names to the output file
with open(snakemake.output.sample_param, 'w') as f:
    f.write(sample_name_str)
