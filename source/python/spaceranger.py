#!/usr/bin/env python3

## Dynamically create the spaceranger count command and then execute command

## Import modules
import subprocess

## Log output and stderr
log_file = open(snakemake.log[0], "w")
sys.stdout = log_file
sys.stderr = log_file

## Handle cases when multiple flow cells were used for sequencing
with open(snakemake.input.sample_param, "r") as file:
    sample_list = file.read().strip()  ## Read and remove any surrounding whitespace
    multi_sample_names = "," in sample_list

# Build the `spaceranger count` command
cmd = [snakemake.params.spaceranger_path, "count",
       f"--id={snakemake.params.run_id}",
       f"--output-dir={snakemake.params.output_dir}",
       f"--transcriptome={snakemake.input.transcriptome_dir}",
       f"--fastqs={snakemake.input.fastqs_dir}",
       f"--slide={snakemake.params.slide}",
       f"--area={snakemake.params.area}",
       f"--probe-set={snakemake.input.probe_set}",
       f"--loupe-alignment={snakemake.input.loupe_alignment}",
       f"--filter-probes=true",
       f"--create-bam=true",
       f"--localcores={snakemake.threads}",
       f"--localmem={snakemake.params.localmem}",
       f"--localvmem={snakemake.params.localvmem}"
      ] 

## Conditionally add --cytaimage if provided
if snakemake.params.cytaimage != "NA":
    cmd.append(f"--cytaimage={snakemake.params.cytaimage}")

## Conditionally add --image if provided
if snakemake.params.image != "NA":
    cmd.append(f"--image={snakemake.params.image}")

## Conditionally add --darkimage if provided
if snakemake.params.darkimage != "NA":
    cmd.append(f"--darkimage={snakemake.params.darkimage}")

## Conditionally add --colorizedimage if provided
if snakemake.params.colorizedimage != "NA":
    cmd.append(f"--colorizedimage={snakemake.params.colorizedimage}")

## Conditionally add --dapi-index if provided
## This will only be provided if darkimage or colorizedimage is present
if snakemake.params.dapi_index != "NA":
    cmd.append(f"--dapi-index={snakemake.params.dapi_index}")

## Add --sample if multiple sample names are present in the FASTQ directory
if multi_sample_names == True:
    cmd.append(f"--sample={sample_list}")
    print("Multiple sample names were identified. Parameter --sample was added to the command.")
else:
    print("Only one sample name was identified. Parameter --sample was not added to the command.")

print(" ".join(cmd))

## Run Space Ranger
subprocess.run(cmd)

## Create an empty file to indicate the Space Ranger run finished without issue
with open(snakemake.output.txt, "w") as file:
    pass
