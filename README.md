## Introduction

kDNA annotation is a python 3.7 package to annotate minicircles of mitochondrial DNA, also known as kinetoplast DNA (kDNA). It identifies cassettes, defined by inverted repeat motifs, and guide RNAs (gRNAs) genes. Alignment of gRNAs with theit cognate mRNAs are found. If gRNA transcriptomics data is available, the expression status of gRNAs are predicted.

## Installation 

It is recommended to have anaconda3 installed. kDNA annotation uses Python3.7 packages pandas, numpy, biopython and yaml. kDNA annotation also requires the C libraries libyaml and the openMP library libgomp1. The package comes with binaries, but if you want to compile the C routines you will need gcc and libyaml-dev.

To clone the latest version of kDNA annotation use

    git clone git@github.com:nicksavill/kDNA-annotation.git

Set up the PATH and PYTHONPATH environment variables in `.bash_profile`: 

    export PYTHONPATH=$HOME/kDNA_annotation:$PYTHONPATH
    export PATH=$HOME/kDNA_annotation/kDNA_annotation:$PATH


## Usage

The complexity and variablility of minicircles within and between species means that minicircle annotation of cassettes and gRNAs cannot be completely automated. User intervention is required throughout the annotation pipeline.

The pipeline commands are contained in the file `share/pipeline.py`. The pipeline is directed by the contents of the configuration file `config.yaml`, an example of which is in the `share/` directory.

kDNA annotation requires, at minimum, the following fasta files:
- minicircle sequences
- unedited mRNA sequences
- edited mRNA sequences

If transcriptomics of gRNAs are available, a gzipped SAM file of transcripts mapped to minicircles is needed.

#### Step 1

For each project set up the following directory structure:
```
    Project_name/
        pipeline.sh
        config.yaml
        In_files/
            minicircles file.fasta
            edited mRNA file.fasta
            unedited mRNA file.fasta
            transcripts file.sam.gz
```
The files `pipeline.sh` and `config.yaml` should be copied from `kDNA_annotation/share/`. In the pipeline two additional directories will be created `Work_files/` and `Annotation/`.

#### Step 2

Edit the `project` parameter in `config.yaml` to the full path name of the project.

#### Step 3

Run `clean_mini_and_maxicircles()` in `pipeline.sh` by uncommenting the call to `clean_mini_and_maxicircles()` and commenting out calls to all the other functions.

For kDNA annotation to work correctly all minicircles should be unique, aligned on conserved sequence block 1 (CSB-1) and have CSB-3 about 100nt downstream from CSB-1. Minicircles should also have the standard name "mO_XXX" where XXX runs from 001 upwards. 

If CSB-1 is missing from any minicircles `clean_mini_and_maxicircles()` will exit with a warning unless `remove minicircles without CSB1` in `config.yaml` is set to yes, then it will automatically remove these minicircles.

Minicircles coataining CSB-1 but not aligned on it will be automatically re-aligned and a warning
