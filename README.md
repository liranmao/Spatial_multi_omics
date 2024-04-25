# Multiplexed spatial mapping of chromatin features, transcriptome, and proteins in tissues

## Introduction
This repository aims to share the raw data processing and visualization codes used in multiplexed spatial mapping of chromatin features, transcriptome, and proteins in tissues paper.


![github_intro](https://github.com/liranmao/Spatial_multi_omics/assets/78578236/2986f184-04c3-4fc2-8166-9e088c83a7cb)


## Data analysis
### Preprocessing
 Next Generation Sequencing (NGS) was performed using the Illumina NovaSeq 6000 sequencer (paired-end 150 bp mode). 
 
In the Data_preprocessing folder, directories beginning with Snakemake_* contain the code for preprocessing different modalities. The preprocessing pipeline utilizes raw FASTQ data as input, where Read 1 comprises genomic sequences, and Read 2 contains the spatial barcodes, specifically Barcode A and Barcode B.

**The preprocessing pipeline we developed using Snakemake workflow management system is in the Data_preprocessing folder. After putting the input files in the correct directory, to run the pipeline, use the command:**

    sbatch Snakemake.sh


**Brief descriptions of preprocessing pipeline in Snakefile:**
##### 1. **Directory and File Setup**
- Automates the creation of directories for storing raw and processed data per sample.
- Lists samples dynamically based on the provided raw data directory.

##### 2. **Reads Filtering**
- `filter_primer`: Utilizes `bbduk.sh` to remove primers and specific sequences from the reads.
- `filter_L1` & `filter_L2`: Further filtering steps target and remove specific linker sequences.

##### 3. **Barcode Processing**
- `bc_process`: Extracts and reformat the data. (BC_process.py)
- `R1_rename`: Renames and reorganizes the processed reads for consistency and further processing.

##### 4. **Barcode Tagging and Matching**
- `taggd`: Corrects barcodes and prepares the data for demultiplexing.
- `add_BC`: Integrat barcode information into sequencing reads based on a pre-generated list of barcodes and their matches to specific reads. (add_BC.py)

##### 5. **Adapter Trimming**
- `adapter_trim`: Trims sequencing adapters from the reads using `trim_galore`, preparing them for alignment.

##### 6. **Sequence Alignment**
- `bwa_align`: Maps reads to a reference genome with `bwa mem`, followed by sorting and indexing the alignments using `samtools`.

##### 7. **Fragment File Generation**
- `fragments`: Transforms BAM files into sorted, compressed, and indexed BED files, ready for downstream analysis.




####  Identify useful pixels (pixel on tissue) from microscope image using Python
See the files in Image_preprocess under Data_preprocessing folder.



### Downstream analysis
All downstream analyses were completed with R language. The package used extensively the functions in Seurat v.4.3.0.1, ArchR v1.0.2, ClusterProfiler v4.8.3, Slingshot v2.2.1, FigR v0.1.0. 

**Brief descriptions of analysis scripts:**

ArchR.R: Analysis of ATAC/histone modifications data.

RNA.R: Analysis of RNA data

GO.R: GO enrichment analysis for marker genes.

Multi_omics_integration.r and Multi_omics_integration_prepare.r: Multi-omics analysis.

RCTD: Cell-type annotation and deconvolution.

FigR: Gene regulation analysis

**Functions:**

spatial_data_visualization.R: Visualize spatially resolved data on tissue sections.



