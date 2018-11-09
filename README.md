# scRNA_pipeline
Snakemake pipeline for Single Cell RNA-seq analysis. Pipeline includes preprocessing,and outlier cells filtering, clustering analysis, annotation and visualization. 


## Prerequisites

 * Software packages
 
   Install all the following R packages:
   
   cidr
   
   cluster
   
   clusterProfiler
   
   dplyr
   
   igraph
   
   MASS
   
   Matrix
   
   org.Hs.eg.db
   
   org.Mm.eg.db
   
   pagoda2
   
   RColorBrewer
   
   scales
   
   scater
   
   scClustViz
   
   scran
   
   Seurat
   
   velocyto.R
   
   viridis

 
 * Please install all the reference files:
   
       mm10:
              transcriptome = "path_to_cell_ranger_download/refdata-cellranger-mm10-1.2.0/mm10"
       gtf = "path_to_cell_ranger_download/refdata-cellranger-mm10-1.2.0/genes/genes.gtf"
              vgtf10x = "path_to_cell_ranger_download/RefGenomes/10X_REF/mm10_rmsk.gtf"
             
       hg19:
              transcriptome = "/path_to_cell_ranger_download/refdata-cellranger-1.1.0/hg19"
              gtf = "path_to_cell_ranger_download/refdata-cellranger-mm10-1.1.0/genes/genes.gtf"
       hg38:
              transcriptome = "path_to_cell_ranger_download/refdata-cellranger-GRCh38-1.2.0"
              gtf = "path_to_cell_ranger_download/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf"
              vgtf10x = "path_to_cellranger_download/RefGenomes/10X_REF/hg38_rmsk.gtf"
        

## Software Installation

   Please download and install all software tools acoording to the specific tool installation instruction. 
   Please download and install the snakemake file from this software package.
   The snakefile will need to be altered to point to the actual package installation locations.

 
   cellranger = "path_to_cell_ranger_download/cellranger-2.2.0/cellranger-cs/2.2.0/bin/cellranger"
   
   velocyto = "path_to_velocyto_tools/velocyto/bin/velocyto"



Once the config file is filled in, then submit submit.sh as a job to run snakemake.
``` 
qsub submit.sh
```

## Running the pipeline
More detailed information about how to run the pipeline can be found on the [wiki](https://github.com/abcsFrederick/scRNA_pipeline/wiki/Single-Cell-RNA-Pipeline-Documentation).

## Contact

  CCRSF_IFX@nih.gov

