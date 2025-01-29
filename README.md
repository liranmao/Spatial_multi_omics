# Multiplexed spatial mapping of chromatin features, transcriptome, and proteins in tissues

## Introduction
This repository aims to share the raw data processing and visualization code used in the **"Multiplexed spatial mapping of chromatin features, transcriptome, and proteins in tissues"** paper. It is now published on Nature Methods ![Multiplexed spatial mapping of chromatin features, transcriptome and proteins in tissues
](https://www.nature.com/articles/s41592-024-02576-0)



![github_intro](https://github.com/liranmao/Spatial_multi_omics/assets/78578236/2986f184-04c3-4fc2-8166-9e088c83a7cb)


## Data analysis
### 1. Preprocessing the sequencing data
 Next Generation Sequencing (NGS) was performed using the Illumina NovaSeq sequencer (paired-end 150 bp mode). 
 
In the Data_preprocessing folder, directories beginning with Snakemake_* contain the code for preprocessing different modalities. The preprocessing pipeline utilizes raw FASTQ data as input, where Read 1 comprises genomic sequences, and Read 2 contains the spatial barcodes, specifically Barcode A and Barcode B.

**The preprocessing pipeline we developed using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system is in the Data_preprocessing folder. After putting the input files in the correct directory, to run the pipeline, use the command:**

    sbatch Snakemake.sh


**Brief descriptions of the preprocessing pipeline in Snakefile:**

(1) Directory and File Setup
- Automates the creation of directories for storing raw and processed data per sample.
- List samples dynamically based on the provided raw data directory.

(2) Reads Filtering
- `filter_primer`: Utilize [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) to filter sequences with primers from the reads.
- `filter_L1` & `filter_L2`: Further filtering steps target specific linker sequences.

(3) Barcode Processing
- `bc_process`: Extract and reformat the data. (BC_process.py)
- `R1_rename`: Rename and reorganize the processed reads for consistency and further processing.

(4) Barcode Tagging and Matching
- `taggd`: Correct barcodes and prepare the data for demultiplexing.
- `add_BC`: Integrate barcode information into sequencing reads based on a pre-generated list of barcodes and their matches to specific reads. (add_BC.py)

(5) Adapter Trimming
- `adapter_trim`: Trims sequencing adapters from the reads using [Trim Galore](https://github.com/FelixKrueger/TrimGalore), preparing them for alignment.

(6) Sequence Alignment
- `bwa_align`: Map reads to a reference genome with [BWA](https://github.com/lh3/bwa), followed by sorting and indexing the alignments using [samtools](https://www.htslib.org/).

(7) Fragment File Generation
- `fragments`: Transforms BAM files into sorted, compressed, and indexed BED files using [sinto toolkit](https://timoast.github.io/sinto/), for downstream analysis. 




###  2. Preprocessing the microscope image

**Identify pixels on tissue from microscope image using Python**

See the files in Image_preprocess under Data_preprocessing folder.



### 3. Downstream analysis
All downstream analyses were completed with R language. The package used extensively the functions in Seurat v.4.3.0.1, ArchR v1.0.2, ClusterProfiler v4.8.3, Slingshot v2.2.1, FigR v0.1.0. To ensure reproducibility, we have provided a comprehensive tutorial demonstrating the joint analysis of multiple modalities and data visualization. You can access the tutorial here: [SpatialMuxSeq vignette](https://rpubs.com/LiranM/SpatialMuxSeq) or in Fig1_joint_analysis.Rmd. The data needed for it can be found on figshare(DOI: [10.6084/m9.figshare.27265410](https://doi.org/10.6084/m9.figshare.27265410.v1)). 

**Brief descriptions of analysis scripts:**

ArchR.R: Analysis of ATAC/histone modifications data.

RNA.R: Analysis of RNA data

GO.R: GO enrichment analysis for marker genes.

Multi_omics_integration.r and Multi_omics_integration_prepare.r: Multi-omics analysis.

RCTD: Cell-type annotation and deconvolution.

FigR: Gene regulation analysis

**Functions:**

spatial_data_visualization.R: Visualize spatially resolved data on tissue sections.


## References

Pengfei Guo, Liran Mao, Yufan Chen, Chin Nien Lee, Angelysia Cardilla, Mingyao Li, Marek Bartosovic, and Yanxiang Deng. "Multiplexed spatial mapping of chromatin features, transcriptome, and proteins in tissues."



