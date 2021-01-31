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

The known CSB-1, 2 and 3 sequences as regular expressions are given in the config file. 

If there are no matches for CSB-1 in a minicircle's sequence this could because the minicircle has not been assembled completely, the regular expression is wrong or that a minor variant of the CSB-1 sequence exists. The program will output possible new variants that it finds at the start the minicircle sequences and then terminate. These variants can be added to the CSB-1 regular expression in the config file. On the other hand, if you want to automatically remove any minicircles without CSB-1 then set `remove minicircles without CSB1` in `config.yaml` to yes.

Minicircles coataining CSB-1 but not aligned on it will be automatically re-aligned and a warning issued.

If CSB-3 does not exist in the correct position relative to CSB-1 a warning is given and the program terminates.

If any minicircle has an invalid name, all minicircles are re-named in the standard format mO_XXX. A log file is written out with the the old and new names of the minicircles. The name of the log file is given by `clean log text file` in the config file.

The maxicircle is also re-named `Maxicircle`.

The cleaned mini and maxicircle fasta files are re-written as `minicircle clean fasta file` and `maxicircle clean fasta file` respectively in the `working directory` as given in the config file.

#### Step 4

Process the edited and unedited mRNA sequences. 

Make sure to comment out `check_minicircles()` and uncomment `mRNA_process()`.

The function `mRNA_process.py` will attempt to clean the inputted edited and unedited sequences given as fasta files in  `in directory`. It is important that the T-stripped edited and unedited sequences are identical otherwise it is usually impossible to reconstruct the insertions and deletions correctly. If this problem occurs, the offending sequences should be corrected. If the mis-match between the sequences occurs in a non-edited region this won't be a problem for gRNA identification. 

Several fasta files are outputted in `working directory`, including positions where deletions occur (`deletions mRNA text file`) and positions where insertions occur which are recorded as a lowercase "u" (`edited_mRNA_small_u.fasta`) or "y" (`edited_mRNA_small_t.fasta`).

#### Step 5

Align the maxicircle (`maxi_align`) and all the minicircles (`mini_align`) to all edited mRNA sequences. This code is written in parellelised C for speed. Both sense and anti-sense strands of the circles are aligned. 

The source code for these binaries are in the `src` directory of the package. If the source code is changed the binaries need re-compiling with the commands

    gcc -O3 -o ../kDNA_annotation/align_maxi align_maxi.c -fopenmp -lyaml
    gcc -O3 -o ../kDNA_annotation/align_mini align_mini.c -fopenmp -lyaml

These programs find all alignments at least as long as `minimum gRNA length` (defined in the config file). As well as finding all gRNAs, these programs also produce many false positives, duplicates and overlapping alignments. In addition, no checks are made for an anhor, nor whether an alignment contains insertions or deletions in the mRNA. These checks are all done at a later step. These programs output to two files `maxi_alignments.txt` and `mini_alignments.txt`. 

#### Step 6

To help in identifying inverted repeats and gRNA genes we first find the high quality gRNA genes; those that we can be quite sure are actual gRNAs rather than false positives. 

High quality gRNAs should be long (>= 40 nt) with a minimal number of mis-matches (<= 1), have a long anchor (>= 8 nt) of only Watson-Crick basepairs. These parameters are defined in the config file.

On running the function `hq_gRNAs.py()` in the pipeline the alignments from Step 5 will be loaded and the high quality gRNAs identified. The number of high quality gRNAs will be printed and the number of minicircles that contain these gRNAs. The gRNA lengths are plotted as a histogram. The histogram should be right-skewed with lengths into the 50s. Also plotted as a histogram are the 5' positions of the high quality gRNAs on the minicircles. These should be clustered into one or more peaks which represent the cassettes. 

If no, or very few high quality gRNAs are found then the filtering parameters in the config file should be weakened and  `hq_gRNAs.py()` re-run.