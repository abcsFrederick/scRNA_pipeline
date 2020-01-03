# iSCAP - An Integrated Single Cell Analysis Pipeline
Snakemake pipeline for Single Cell RNA-seq analysis. The integrated Single Cell Analysis Pipeline (iSCAP) enables combining different methods including quality control, preprocessing, normalization, batch effect correction, clustering, marker gene and subpopulation identification, cell type annotation, multi-sample differential expression analysis. 

The iSCAP allows users to select and apply their preferred combination of the software tools supported in the pipeline and define sample groups for comparison analysis.  iSCAP also includes gene set enrichment and network analysis, it generates HTML-based nozzle reports which allow users to browse the results obtained in each stage of the pipeline. 


## Prerequisites

 * Software packages
 
   Install all the following R packages:
   
   ```
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
   reshape2
   scales
   scater
   scClustViz
   scran
   Seurat (version 2)
   SingleR
   URD
   velocyto.R
   viridis
   ```
 
 * Please install all the reference files:
   
       mm10:
              transcriptome = "path_to_cell_ranger_download/refdata-cellranger-mm10-3.0.0"
       gtf = "path_to_cell_ranger_download/refdata-cellranger-mm10-3.0.0/genes/genes.gtf"
              vgtf10x = "path_to_cell_ranger_download/RefGenomes/10X_REF/mm10_rmsk.gtf"
             
       hg19:
              transcriptome = "path_to_cell_ranger_download/refdata-cellranger-hg19-3.0.0/"
              gtf = "path_to_cell_ranger_download/refdata-cellranger-hg19-3.0.0/genes/genes.gtf"
       hg38:
              transcriptome = "path_to_cell_ranger_download/refdata-cellranger-GRCh38-3.0.0"
              gtf = "path_to_cell_ranger_download/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
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

