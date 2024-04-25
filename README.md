# Multiplexed spatial mapping of chromatin features, transcriptome, and proteins in tissues

## Introduction
This repository aims to share the raw data processing and visualization codes used in multiplexed spatial mapping of chromatin features, transcriptome, and proteins in tissues paper.


![github_intro](https://github.com/liranmao/Spatial_multi_omics/assets/78578236/2986f184-04c3-4fc2-8166-9e088c83a7cb)


## Data analysis
### Preprocessing
 Next Generation Sequencing (NGS) was performed using the Illumina NovaSeq 6000 sequencer (paired-end 150 bp mode). In Data_preprocessing folder, the folders with the name starting with Snakemake_* is the code for preprocessing different modality. The input for the preprocessing pipeline is raw fastq data with read 1 contains the genome sequences and read 2 contains the spatial Barcode A and Barcode B. 

**The preprocessing pipeline we developed using Snakemake workflow management system is in the Data_preprocessing folder. After putting the input files in the correct directory, to run the pipeline, use the command:**

    sbatch Snakemake.sh

**Brief descriptions of preprocessing pipeline:**
###### 1. **Directory and File Setup**
- Automates the creation of directories for storing raw and processed data per sample.
- Lists samples dynamically based on the provided raw data directory.

#### 2. **Reads Filtering**
- `filter_primer`: Utilizes `bbduk.sh` to remove primers and specific sequences from the reads.
- `filter_L1` & `filter_L2`: Further filtering steps target and remove specific linker sequences.

#### 3. **Barcode Processing**
- `bc_process`: Extracts and processes barcode information from the reads for spatial identification.
- `R1_rename`: Renames and reorganizes the processed reads for consistency and further processing.

#### 4. **Barcode Tagging and Matching**
- `taggd`: Corrects barcodes and prepares the data for demultiplexing.
- `add_BC`: Attaches barcode tags back to the reads, facilitating tracking and analysis.

#### 5. **Adapter Trimming**
- `adapter_trim`: Trims sequencing adapters from the reads using `trim_galore`, preparing them for alignment.

#### 6. **Sequence Alignment**
- `bwa_align`: Maps reads to a reference genome with `bwa mem`, followed by sorting and indexing the alignments using `samtools`.

#### 7. **Fragment File Generation**
- `fragments`: Transforms BAM files into sorted, compressed, and indexed BED files, ready for downstream analysis.




    

#### Spatial_ATAC-seq
##### 1.Raw Fastq data

##### 2. Reformat raw Fastq file to Cell Ranger ATAC format (10x Genomics)
**Raw read 1 -> New Read 1 + New Read 2**
- New Read 1: contains the genome sequences
- New Read 2: contains the spatial Barcode A and Barcode B

**Raw read 2 -> New Read 3**

Reformatting raw data was implemented by BC_process.py in the Data_preprocessing folder.


##### 3. Sequence alignment and generation of fragments file
The reformated data was processed using Cell Ranger ATAC v1.2 with following references:
Mouse reference (mm10):

    curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-mm10-1.2.0.tar.gz

Human reference (GRCh38):

    curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-1.2.0.tar.gz

**A preprocessing pipeline we developed using Snakemake workflow management system is in the Data_preprocessing folder. To run the pipeline, use the command:**

    sbatch Snakemake.sh

#### Spatial_RNA-seq

**1. Raw Fastq data processing using ST pipeline and generate expression matrix**

The DBiT-seq Raw fastq file
Read 1: Contains the cDNA sequence
Read 2: Contains the spatial Barcode A, Barcode B and UMIs

**2.Reformat Fastq Read 2 file**
To reformat the Raw data, run the fastq_process.py in Rawdata_processing folder and gzip the resulted fastq file to save space:

    python fastq_process.py
    gzip sample_R2_processed.fastq

The reformated data was processed following ST pipeline.  

**3、Run ST pipeline**
Run st_pipeline.sh to start the ST pipeline: The input is processed_R2.fastq.gz and Raw R1.fastq.gz. It also requires a "spatial_barcodes_index.txt" to decode the spatial location information. Genome references and annotatation files were aslo needed.

    #!/bin/bash
    # FASTQ reads
    FW=$tmp/${sample}_R2_processed.fastq
    RV=$tmp/${sample}_raw_qc_R1.fastq.gz
    # References for mapping, annotation and nonRNA-filtering
    MAP=/Dropseq_Alignment_References/mm10/
    ANN=/Dropseq_Alignment_References/mm10/mm10.gtf 
    CONT=/Spatial_omics_references/mouse/GRCm38_86/ncRNA/StarIndex/

    # Barcodes settings
    ID=/useful/spatial_barcodes.txt

    # Output folder and experiment name
    OUTPUT=../output/
    mkdir -p $OUTPUT

    TMP_ST=$OUTPUT/tmp
    mkdir -p $TMP_ST

    # Running the pipeline
    st_pipeline_run.py \
      --output-folder $OUTPUT \
      --ids $ID \
      --ref-map $MAP \
      --ref-annotation $ANN \
      --expName $sample \
      --htseq-no-ambiguous \
      --verbose \
      --log-file $OUTPUT/${sample}_log.txt \
      --demultiplexing-kmer 5 \
      --threads 20 \
      --temp-folder $TMP_ST \
      --no-clean-up \
      --umi-start-position 16 \
      --umi-end-position 26 \
      --demultiplexing-overhang 0 \
      --min-length-qual-trimming 20 \
      $FW $RV

**4.Convert Ensemble to Gene Names**
Then, Run converttoname.sh to annotate the resulting sample_stdata.tsv.
    
    #!/bin/bash
    tsv_E=$OUTPUT/${sample}_stdata.tsv
    path_to_annotation_file=/Dropseq_Alignment_References/mm10/mm10.gtf

    convertEnsemblToNames.py --annotation $path_to_annotation_file --output $OUTPUT/${sample}_stdata_names.tsv $tsv_E

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



