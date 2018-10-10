# scRNA_pipeline
Single Cell RNA-seq Analysis Pipeline


Copy or generate a config.py file in the folder that snakemake will be run in.
The config.py file is composed of the following lines:
 
unaligned="demultiplex_path"
analysis="analysis_folder"
ref="reference_genome"
numcells="3000,3000"
 
Unaligned - the direct path to the folder containing the fastq files for each sample
Analysis - the location that the analysis will take place in
Reference genome - the reference genome (mm10 or hg38)
NumCells - the estimated number of cells for each sample, with the default being 3000

Once the config file is filled in, then submit submit.sh as a job to run snakemake.
 
qsub submit.sh
