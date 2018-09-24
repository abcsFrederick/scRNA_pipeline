#Snakefile for single cell pipeline 
#Ashley Walton
#CCR-SF IFX 

import config
from snakemake.utils import R
import glob
import os.path 
import os


genomename = ""
if config.ref == "mm10":
	transcriptome = "/installed_tools/GemCode/cell_ranger_download/refdata-cellranger-mm10-1.2.0"
	gtf = "/installed_tools/GemCode/cell_ranger_download/refdata-cellranger-mm10-1.2.0/genes/genes.gtf"
	vgtf10x = "/installed_tools/RefGenomes/10X_REF/mm10_rmsk.gtf"
	genomename = "mm10"
if config.ref == "hg19":
	transcriptome = "/installed_tools/RefGenomes/10X_REF/refdata-cellranger-1.1.0/hg19"
	genomename = "hg19"
if config.ref == "hg38":
	transcriptome = "/installed_tools/GemCode/cell_ranger_download/refdata-cellranger-GRCh38-1.2.0"
	gtf = "/installed_tools/GemCode/cell_ranger_download/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf"
	vgtf10x = "/is2/projects/CCR-SF/active/RefGenomes/10X_REF/hg38_rmsk.gtf"
	genomename = "GRCh38"

rscripts = "/installed_tools/scripts/R_scripts"
bclpathbase = "/installed_tools/illumina/RawData_NextSeq/"
fastqpath = config.unaligned
analysis = config.analysis

cellranger = "/installed_tools/GemCode/cell_ranger_download/cellranger-2.2.0/cellranger-cs/2.2.0/bin/cellranger"
velocyto = "/installed_tools/velocyto/bin/velocyto"

numcell = config.numcells
numcells = config.numcells.split(',')
sample = [os.path.basename(file).split('.')[0] for file in glob.glob(fastqpath+'/*')]
samps = []
for item in sample:
        newvar = item.split("_R1")
        othervar = item.split("_R2")
        samps.append(newvar[0])
new = []
for item in samps:
        if '_R2_' not in item:
                new.append(item)
samples = [s.replace('Sample_', '') for s in new]
samples = sorted(samples)

dict2 = dict(zip(samples,numcells))

rule all:
	input: expand("run_{sample}_10x_cellranger_count.err", sample=samples), expand("{sample}/seurat/minRes.txt", sample=samples), expand("{sample}/clusterpro/groupGO.pdf", sample=samples), expand("{sample}/velocyto/{sample}.loom", sample=samples), expand("{sample}/velocyto/velocyto.pdf", sample=samples), expand("{sample}/scran/CellCycle.pdf", sample=samples)

rule count: 
	output: err = "run_{sample}_10x_cellranger_count.err", log ="run_{sample}_10x_cellranger_count.log"
	params: batch = "-l nodes=1:ppn=16,mem=96gb", prefix = "{sample}", prefix2 = fastqpath + "/{sample}", numcells = lambda wildcards:dict2[wildcards.sample]
	shell: "module load bcl2fastq2; {cellranger} count --id={params.prefix} --fastqs={params.prefix2} --expect-cells={params.numcells} --transcriptome={transcriptome} 2>{output.err} 1>{output.log}" 

rule seurat:
	input: "run_{sample}_10x_cellranger_count.err"
	output: rds = "{sample}/seurat/seur_10x_cluster_object.rds", mr = "{sample}/seurat/minRes.txt", gm = "{sample}/seurat/Filtered_Gene_Expression_Matrix.csv"
	log: "{sample}/seurat/seurat.log"
	params: batch = "-l nodes=1:ppn=8,mem=48gb", prefix = "{sample}", prefix2 = "{sample}/outs/filtered_gene_bc_matrices/"+genomename, outdir = "{sample}/seurat"
	shell: "set -xv; export R_LIBS=/opt/nasapps/applibs/r-3.5.0_libs/;/opt/nasapps/development/R/3.5.0/bin/R  --no-save --args {params.outdir} {params.prefix2} {params.prefix} {config.ref} <{rscripts}/Final_Seurat.R >& {log}"
	
rule clusterpro:
	input: "{sample}/seurat/minRes.txt"
	output: pdf = "{sample}/clusterpro/groupGO.pdf"
	log: "{sample}/clusterpro/clusterpro.log"
	params: batch = "-l nodes=1:ppn=8", outdir = "{sample}/clusterpro", prefix = "{sample}", prefix2 = "{sample}/seurat"
	shell: "set -xv; export R_LIBS=/opt/nasapps/applibs/r-3.5.0_libs/;/opt/nasapps/development/R/3.5.0/bin/R --no-save --args {params.outdir} {config.ref} {params.prefix2} {input} <{rscripts}/clusterProfiler_seurat.R >& {log}"	

rule pvelocyto:
    input: "run_{sample}_10x_cellranger_count.err"
    output: out = "{sample}/velocyto/{sample}.loom"
    log: "{sample}/velocyto/pvelocyto.log"
    params: batch = "-l nodes=1:ppn=16", prefix = "{sample}"
    shell: "export PATH=/installed_tools/Anaconda/3.6/install/bin:$PATH; export PATH=/opt/nasapps/development/samtools/1.7/bin/:$PATH; export PYTHONPATH=/installed_tools/velocyto/lib/python3.6/site-packages/; {velocyto} run10x -@ 16 --samtools-memory 48000 -m {vgtf10x} {params.prefix} {gtf} >& {log}"	
	
rule svelocyto:
    input: loom = "{sample}/velocyto/{sample}.loom", rds = "{sample}/seurat/seur_10x_cluster_object.rds"
    output: "{sample}/velocyto/velocyto.pdf" 
    log: "{sample}/velocyto/svelocyto.log"
    params: batch = "-l nodes=1:ppn=8,mem=64gb", prefix = "{sample}/velocyto"
    shell: "export R_LIBS=/opt/nasapps/applibs/r-3.5.0_libs/;/opt/nasapps/development/R/3.5.0/bin/R --no-save --args {params.prefix} {input.loom} {input.rds} < {rscripts}/velocyto_seurat.R >& {log}"
	
rule scran:
    input: "{sample}/seurat/Filtered_Gene_Expression_Matrix.csv"
    output: cellcycle = "{sample}/scran/CellCycle.pdf", lastplot = "{sample}/scran/ExpressionPlot.pdf"
    log: "{sample}/scran/scranp.log"
    params: batch = "-l nodes=1:ppn=8,mem=64gb", outdir = "{sample}/scran"
    shell: "set -xv; export R_LIBS=/opt/nasapps/applibs/r-3.5.0_libs; /opt/nasapps/development/R/3.5.0/bin/R --no-save --args {input} {params.outdir} {config.ref} <{rscripts}/scran_cellcyclePlot.R >& {log}"
	
onsuccess:
    success = "Yes"
    al = analysis
    shell(r'''python /is2/projects/CCR-SF/active/Software/scripts/bin/sendmail.py {success} {al}''')

onerror:
    success = "No"
    al = analysis
    shell(r'''python /is2/projects/CCR-SF/active/Software/scripts/bin/sendmail.py {success} {al}''')