#############################################################################################################
# pipeline_template
#
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
#
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#
# Program:  Snakefile
# Version:  0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  Main snakefile for workflow template.
# Input:    Sample fastq files stored in a cloud location (google cloud, aws)
# Output:   'sample_metadata.csv','rseqc/bam_qc_parsed.tab', 'rseqc/bam_gc_parsed.tab'
#############################################################################################################

## Import modules
import shutil
import logging as _logging
import psutil
import os
from pathlib import Path
from box import Box
import yaml
import pandas as pd
from typing import Literal
from typing import Optional

wildcard_constraints:
    sample="[^_]+"


##############################
#       CONFIGURATION        #
##############################
# Specify YAML config file location
configfile: "config/config.yaml"

# Directories
# working output dir
PREDIR     = config["predir"]
# source dir for supporting scripts
SOURCEDIR  = config["srcdir"]
# analysis data results and reporting within working dir
DATADIR    = PREDIR+'/analysis/data'
REPDIR     = PREDIR+'/analysis/report'
# path location for the sample_metadata.csv
METADATA_CSV = os.path.join(os.path.dirname(SOURCEDIR), config["sample_metadata"])

# Use the source dir to import helper modules
try:
    sys.path.append(SOURCEDIR+'/python')
    import trimadapters
    import getfile
    import putfile
    import utils
    import time
except:
    print("The srcdir value in the config file has not been properly configured. \
           Please configure the config/config.yaml file and try again.")

#added back in for to_log and to_benchmark functions
include: "./rules/common.smk"


## create file accessor
paths = create_path_accessor()

## read in reference genome locations file
reference_df = pd.read_table(config["reference"], sep=",")

## TODO: Check to see if sample_metadata_df can be completely replaced by sample_df
## read in sample metadata file
sample_metadata_df = pd.read_table(config["sample_metadata"], sep=",", keep_default_na=False, comment = "#")

## Add information on Visium type to the sample_metadata_df
## Refer to https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/image-slide-parameter
## Slide starts with "V" = Visium Spatial Gene Expression or Visium CytAssist Spatial Gene Expression
## Slide starts with "H" = Visium HD
type_map = {"V": "sge", "H": "hd"}
## Apply mapping and force an error for invalid cases for any NaN values
sample_metadata_df["visium_type"] = sample_metadata_df["slide"].str[0].map(type_map)
if sample_metadata_df["visium_type"].isna().any():
    invalid_slides = sample_metadata_df.loc[sample_metadata_df["visium_type"].isna(), "slide"].tolist()
    raise ValueError(f"Invalid slide values: {invalid_slides}. Expected to start with 'V' or 'H'.")
print(sample_metadata_df)

## Extract the file extension for 'image' and 'cytaimage' and information as additional columns
#sample_df = sample_metadata_df.set_index('samid')
#sample_df['image_ext'] = sample_df['image'].apply(lambda x: os.path.splitext(x)[1] if isinstance(x, str) else '')
#sample_df['cytaimage_ext'] = sample_df['cytaimage'].apply(lambda x: os.path.splitext(x)[1] if isinstance(x, str) else '')
#print(sample_df)
print(sample_metadata_df)

## Reference genome files for Space Ranger
REF_GENOME = reference_df.loc[reference_df["ref_file_name"]=="ref_genome", "google_bucket_URI"].item()

## Probe sets for Space Ranger
PROBE_V1 = reference_df.loc[reference_df["ref_file_name"]=="probe_v1", "google_bucket_URI"].item()
PROBE_V2 = reference_df.loc[reference_df["ref_file_name"]=="probe_v2", "google_bucket_URI"].item()

## Space Ranger tar file
SR_TAR = reference_df.loc[reference_df["ref_file_name"]=="sr_tar", "google_bucket_URI"].item()

## Return the input for rules (tabula_sapiens, custom, clustering) depending on which optional
## rules were performed previously. Optional rules in order of execution are
## cell_cycle_scoring, azimuth, tabula_sapiens, and custom.
def rule_input(cell_cycle_score: Literal["TRUE", "FALSE"],
               azimuth: Literal["TRUE", "FALSE"],
               ts: Optional[Literal["TRUE", "FALSE"]],
               custom: Optional[Literal["TRUE", "FALSE"]], paths: Box) -> Path:
    ## If custom is performed, its RDS is returned since it is run last
    if custom == "TRUE":
        return paths.cell_type_annotation.rds_custom
    ## If tabula_sapiens is performed, and custom is NOT performed
    if ts == "TRUE":
        return paths.cell_type_annotation.rds_ts
    ## If azimuth is performed, and tabula_sapiens and custom are NOT performed
    if azimuth == "TRUE":
        return paths.cell_type_annotation.rds_azimuth
    ## If cell_cycle_score is performed, and azimuth, tabula_sapiens, and custom are NOT performed
    if cell_cycle_score == "TRUE":
        return paths.cell_cycle_scoring.rds
    ## Default case if none of the rules were performed (cell_cycle_score, azimuth, tabula_sapiens, or custom)
    return paths.dimensionality_reduction.rds

## Define input RDS for rule tabula_sapiens
TS_INPUT = rule_input(config["cell_cycle_scoring"]["cell_cycle_score"],
                      config["cell_type_annotation"]["azimuth"], None, None, paths)

## Define input RDS for rule custom
CUSTOM_INPUT = rule_input(config["cell_cycle_scoring"]["cell_cycle_score"],
                          config["cell_type_annotation"]["azimuth"],
                          config["cell_type_annotation"]["ts"], None, paths)

## Define input RDS for rule clustering
CLUSTERING_INPUT = rule_input(config["cell_cycle_scoring"]["cell_cycle_score"],
                              config["cell_type_annotation"]["azimuth"],
                              config["cell_type_annotation"]["ts"],
                              config["cell_type_annotation"]["custom"], paths)

print("Input for rule tabula_sapiens: " + TS_INPUT)
print("Input for rule custom: " + CUSTOM_INPUT)
print("Input for rule clustering: " + CLUSTERING_INPUT)

## Generate the output path for the Tabula Sapiens reference file
def ref_path(gcp_uri: str) -> Path:
    name = os.path.basename(gcp_uri)
    path = Path(PREDIR) / "ref_files" / name
    return path

## Generate the output path for the reference file URI. The reference dataframe is not used to hold these references
## since they can vary and are not set.
TS_PATH = ref_path(config["cell_type_annotation"]["ts_ref"])
CUSTOM_PATH = ref_path(config["cell_type_annotation"]["custom_ref"])
CUSTOM_COUNTS = ref_path(config["cell_type_annotation"]["custom_counts"])
CUSTOM_METADATA = ref_path(config["cell_type_annotation"]["custom_metadata"])

print(TS_PATH)

# Sample info
## List of samples to process
SAMID = utils.toList(sample_metadata_df['samid'])



## List of sample(s) that have a cytaimage
CYTAIMAGE = sample_metadata_df.loc[
        sample_metadata_df["cytaimage"].notna() & sample_metadata_df["cytaimage"].str.strip().ne(""),
        "samid"].tolist()
#print(f"Cytaimage: {CYTAIMAGE}")

## List of sample(s) that have an image
IMAGE = sample_metadata_df.loc[
        sample_metadata_df["image"].notna() & sample_metadata_df["image"].str.strip().ne(""),
        "samid"].tolist()
#print(f"Image: {IMAGE}")

## List of sample(s) that have a darkimage
DARKIMAGE = sample_metadata_df.loc[
        sample_metadata_df["darkimage"].notna() & sample_metadata_df["darkimage"].str.strip().ne(""),
        "samid"].tolist()
#print(f"Darkimage:  {DARKIMAGE}")

## List of sample(s) that have a colorizedimage
COLORIZEDIMAGE = sample_metadata_df.loc[
        sample_metadata_df["colorizedimage"].notna() & sample_metadata_df["colorizedimage"].str.strip().ne(""),
                "samid"].tolist()
#print(f"Colorizedimage: {COLORIZEDIMAGE}")

# Set workflow working (output) dir
workdir: PREDIR



##############################
#       SET GLOBAL VARS      #
##############################
# Workflow info
## number of cores dedicated to run
NCORES = int(config["ncores"])
## initial sub folders
SUBDIRS = 'benchmark log info progress ref_files input analysis analysis/data analysis/report'

# Cloud options retrieving files and archiving results
CLOUD = config["cloud_prog"]
ARCHIVE = config["archive_bucket"]
DOARCHIVE = len(ARCHIVE) != 0

# Set up logging to log file
LOG_FILE  = config["log_file"]
_logging.basicConfig(level=_logging.INFO,
                    format='[pipeline_template] %(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    handlers=[_logging.FileHandler(LOG_FILE)])

## TODO: Rename "image_path" since the function is used for more than just retrieving path of the images.
## Generate image PATHs for Space Ranger image and cytaimage parameters. The file types can vary.
## Therefore, the PATH needs to be created dynamically based on what the file
## type is.
def image_path(sample, image_param):
    #print(sample)
    #print(image_param)
    image_uri = sample_metadata_df.loc[sample_metadata_df['samid'] == sample, image_param].values[0]

    if image_param in {"cytaimage", "image", "darkimage", "colorizedimage"}:
        image_extension = Path(image_uri).suffix
        image_path = f"{PREDIR}/input/{sample}_inputs/{sample}_{image_param}{image_extension}"
        #print(image_uri)
        #print(image_path)
        return image_uri, image_path
    else:
        return image_uri

#value = image_path("breast.sge1", "image")[1]
#print(value)

################################
#     DEFINE TARGET OUTPUT     #
################################
OUTPUT = [paths.ref_files.refdata,
          PREDIR+"/progress/spaceranger_env.done",
          expand(paths.input.fastqs_dir, sample=SAMID),
          expand(paths.input.cytaimage, sample=CYTAIMAGE), ## This needs to always be in the OUTPUT
          expand(paths.input.image, sample=IMAGE), ## This needs to always be in the OUTPUT
          expand(paths.input.darkimage, sample=DARKIMAGE), ## This needs to always be in the OUTPUT
          expand(paths.input.colorizedimage, sample=COLORIZEDIMAGE), ## This needs to always be in the OUTPUT
          expand(paths.input.loupe_json, sample=SAMID),
          expand(paths.input.sample_param, sample=SAMID),
          expand(paths.spaceranger.txt, sample=SAMID),
          expand(paths.qc.rds, sample=SAMID),
          expand(paths.qc.rds_qc, sample=SAMID),
          expand(paths.filter.rds, sample=SAMID),
          expand(paths.normalization.rds, sample=SAMID),
          expand(paths.dimensionality_reduction.rds, sample=SAMID),
          expand(paths.clustering.rds, sample=SAMID),
          expand(paths.spatial_variable_features.rds, sample=SAMID),
          expand(paths.report.html, sample=SAMID),
          #paths.merge.rds
         ]

if config["cell_cycle_scoring"]["cell_cycle_score"] == 'TRUE':
    OUTPUT.append(expand(paths.cell_cycle_scoring.rds, sample=SAMID))

if config["cell_type_annotation"]["azimuth"] == 'TRUE':
    OUTPUT.append(expand(paths.cell_type_annotation.rds_azimuth, sample=SAMID))

if config["cell_type_annotation"]["ts"] == 'TRUE':
    OUTPUT.append(expand(paths.cell_type_annotation.rds_ts, sample=SAMID))

if config["cell_type_annotation"]["custom"] == 'TRUE':
    OUTPUT.append(expand(paths.cell_type_annotation.rds_custom, sample=SAMID))

## Output for the 10X workflow. This will skip the steps up to and including the STARsolo alignment.
## The input is the 10X output folder that contains the barcodes.tsv, features.tsv, and matrix.mtx
#start_index = OUTPUT.index(expand(paths.qc.vln, sample=SAMID))
#OUTPUT_10X = OUTPUT[start_index:]

#print(OUTPUT)
#print(OUTPUT_10X)

#########################################
#    Define any onstart or onsuccess    #
#########################################
onsuccess:
    ## Copy sample_metadata.csv to the PREDIR
    shell('cp '+SOURCEDIR+'/../'+config["sample_metadata"]+' '+PREDIR)

    #if config["workflow"] != "10x":
        ## Merge sample rseqc results into single result files
        #merged_results = utils.mergeRSEQC(SOURCEDIR)

        ## Copy some results to analysis data dir
        #[shutil.copy2(x, DATADIR) for x in merged_results]

        ## knit rmarkdown html report
        #shell('Rscript --vanilla '+SOURCEDIR+'/r/run-report.r '+SOURCEDIR+'/r/cidc_atac-report-slidy.Rmd '+PREDIR+' '+DATADIR+'/../report')
        #merged_results.append('analysis/report/cidc_atac-report-slidy.html')

    ## Upload main results if needed
    #if DOARCHIVE:
        #[putfile.upload(file=x, destination=ARCHIVE, prog=CLOUD) for x in merged_results]

    shell("echo 'Pipeline complete!'")



################################
#   PASS OUTPUT TO all RULE    #
################################
rule all:
    input:
        OUTPUT_10X if config["workflow"] == "10x" else OUTPUT

################################
#        PIPELINE RULES        #
################################
include: "./rules/initialization.smk"
include: "./rules/ingest.smk"
include: "./rules/mapping.smk"
include: "./rules/qc.smk"
include: "./rules/filter.smk"
include: "./rules/normalization.smk"
include: "./rules/dimensionality_reduction.smk"
include: "./rules/cell_cycle_scoring.smk"
include: "./rules/cell_type_annotation.smk"
include: "./rules/clustering.smk"
include: "./rules/spatial_variable_features.smk"
include: "./rules/report.smk"
include: "./rules/merge.smk"
