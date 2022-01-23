## Introduction

**vibbits/rnaseq-editing** is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation followed by a prediction step of editing sites using RDDpred.

The pipeline is largely based on the [nf-core RNAseq pipeline](https://nf-co.re/rnaseq/).

The initial nf-core pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Adapter and quality trimming ([`Trimmomatics`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
4. Use of STAR for multiple alignment and quantification: [`STAR`](https://github.com/alexdobin/STAR)
5. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
6. Prediction of editing sites using RDDpred ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
7. Extensive quality control:
    1. [`RSeQC`](http://rseqc.sourceforge.net/)
    2. [`Qualimap`](http://qualimap.bioinfo.cipf.es/)
    3. [`dupRadar`](https://bioconductor.org/packages/release/bioc/html/dupRadar.html)
8. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) on a Linux operating system.
   Note: This pipeline does not currently support running with macOS.

3. Download the pipeline via git clone, download the associated training data files for RDDpred into the assets folder, download the reference data to 

    ```console
    git clone https://github.com/vibbits/rnaseq-editing.git
    cd $(pwd)/rnaseq-editing/assets
    # download training data file for RDDpred
    wget -c 
    # download reference data for your genome, we provide genome and indexed genome for STAR 2.7.3a
    
    ```

    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis using Docker locally!

    ```console
    nextflow run vibbits/rnaseq-editing \
        --input samplesheet.csv \
        --genome hg19 \
        -profile docker
    ```

    * An executable Python script called [`fastq_dir_to_samplesheet.py`](https://github.com/nf-core/rnaseq/blob/master/bin/fastq_dir_to_samplesheet.py) has been provided if you would like to auto-create an input samplesheet based on a directory containing FastQ files **before** you run the pipeline (requires Python 3 installed locally) e.g.

        ```console
        wget -L https://raw.githubusercontent.com/nf-core/rnaseq/master/bin/fastq_dir_to_samplesheet.py
        ./fastq_dir_to_samplesheet.py <FASTQ_DIR> samplesheet.csv --strandedness reverse
        ```

    * The final analysis has been executed on the Azure platform using Azure Kubernetes Services (AKS). AKS has to be set up on the Azure platform by defining a standard node pool called sys next to the scalable node pool cpumem using Standard_E8ds_v4 as node size for calculation.
      Furthermore, persistent volume claims (PVCs) have been setup for input and work folders of the nextflow runs. In the PVC `input` the reference data as well as the fastqc files have been stored where the PVC `work`, the temporary nextflow files for the individual runs as well as the output files have been stored.
    * The config file for the final execution run for [RNAseq editing for the human samples and reference genome hg19](https://github.com/vibbits/rnaseq-editing/blob/master/nextflow.config.as-executed).    

## Documentation

The nf-core/rnaseq pipeline comes with documentation about the pipeline [usage](https://nf-co.re/rnaseq/usage), [parameters](https://nf-co.re/rnaseq/parameters) and [output](https://nf-co.re/rnaseq/output).

## Credits
These scripts were written to provide a reproducible data analysis pipeline until the downstream processing using dedicated R scripts for exploratory analysis and plotting. The general structure of pipeline is based on the data analysis steps of the our recent paper [ADAR1 interaction with Z-RNA promotes editing of endogenous double-stranded RNA and prevents MDA5-dependent immune activation](https://pubmed.ncbi.nlm.nih.gov/34380029/).

Note: The nf-core scripts this pipeline is based on were originally written for use at the [National Genomics Infrastructure](https://ngisweden.scilifelab.se), part of [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

The RNAseq pipeline was re-written in Nextflow DSL2 by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) from [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

## Citations

The `nf-core` publication is cited here as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
