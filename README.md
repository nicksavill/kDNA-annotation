## Introduction

kDNA annotation is a python 3.7 package to annotate minicircles of mitochondrial DNA, also known as kinetoplast DNA (kDNA). It identifies cassettes, defined by inverted repeat motifs, and guide RNAs (gRNAs) genes. Alignment of gRNAs with their cognate mRNAs are found. If gRNA transcriptomics data is available, the expression status of gRNAs are predicted.

## Installation 

It is recommended to have anaconda3 installed. kDNA annotation uses Python3.7 packages pandas, numpy, biopython and yaml. kDNA annotation also requires the C libraries libyaml and the openMP library libgomp1. The package comes with binaries, but if you want to compile the C routines you will need gcc and libyaml-dev.

To clone the latest version of kDNA annotation use

    git clone git@github.com:nicksavill/kDNA-annotation.git

Set up the PATH and PYTHONPATH environment variables in `.bash_profile`: 

    export PYTHONPATH=$HOME/kDNA_annotation:$PYTHONPATH
    export PATH=$HOME/kDNA_annotation/kDNA_annotation:$PATH

The pipeline also uses `meme` from the https://meme-suite.org/meme/. You can either use their server or download and install the software.

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
        pipeline.py
        config.yaml
        In_files/
            minicircles file.fasta
            edited mRNA file.fasta
            unedited mRNA file.fasta
            transcripts file.sam.gz
```
The files `pipeline.py` and `config.yaml` should be copied from `kDNA_annotation/share/`. In the pipeline two additional directories will be created `Work_files/` and `Annotation/`.

#### Step 2

Edit the `project` parameter in `config.yaml` to the full path name of the project.

#### Step 3

Clean the minicircle and maxicircle fasta files with `clean_mini_and_maxicircles()`.

Run `clean_mini_and_maxicircles()` in `pipeline.py` by un-commenting the call to `clean_mini_and_maxicircles()` and commenting out calls to all the other functions.

For kDNA annotation to work correctly all minicircles should be unique, aligned on conserved sequence block 1 (CSB-1) and have CSB-3 about 100nt downstream from CSB-1. Minicircles should also have the standard name "mO_XXX" where XXX runs from 001 upwards. 

The known CSB-1, 2 and 3 sequences as regular expressions are given in the config file. 

If there are no matches for CSB-1 in a minicircle's sequence this could because the minicircle has not been assembled completely, the regular expression is wrong or that a minor variant of the CSB-1 sequence exists. The program will output possible new variants that it finds at the start the minicircle sequences and then terminate. These variants can be added to the CSB-1 regular expression in the config file. On the other hand, if you want to automatically remove any minicircles without CSB-1 then set `remove minicircles without CSB1` in `config.yaml` to yes.

Minicircles containing CSB-1 but not aligned on it will be automatically re-aligned and a warning issued.

If CSB-3 does not exist in the correct position relative to CSB-1 a warning is given and the program terminates.

If any minicircle has an invalid name, all minicircles are re-named in the standard format mO_XXX. A log file is written out with the the old and new names of the minicircles. The name of the log file is given by `clean log text file` in the config file.

The maxicircle is also re-named `Maxicircle`.

The cleaned mini and maxicircle fasta files are re-written as `minicircle clean fasta file` and `maxicircle clean fasta file` respectively in the `working directory` as given in the config file.

#### Step 4

Process the edited and unedited mRNA sequences with `mRNA_process()`.

Make sure to comment out `clean_mini_and_maxicircles()` and uncomment `mRNA_process()`.

The function `mRNA_process.py` will attempt to clean the inputted edited and unedited sequences given as fasta files in  `in directory` defined in the config file. It is important that the T-stripped edited and unedited sequences are identical otherwise it is usually impossible to reconstruct the insertions and deletions correctly. If this problem occurs, the offending sequences should be corrected. If the mis-match between the sequences occurs in a non-edited region this won't be a problem for gRNA identification. 

Several fasta files are outputted in `working directory`, including positions where deletions occur (`deletions mRNA text file`) and positions where insertions occur which are recorded as a lowercase "u" (`edited_mRNA_small_u.fasta`) or "y" (`edited_mRNA_small_t.fasta`).

#### Step 5

Align the maxicircle (`maxi_align`) and all the minicircles (`mini_align`) to all edited mRNA sequences. This code is written in parallelised C for speed. Both sense and anti-sense strands of the circles are aligned. 

The source code for these binaries are in the `src` directory of the package. If the source code is changed the binaries need re-compiling with the commands

    gcc -O3 -o ../kDNA_annotation/align_maxi align_maxi.c -fopenmp -lyaml
    gcc -O3 -o ../kDNA_annotation/align_mini align_mini.c -fopenmp -lyaml

These programs find all alignments at least as long as `minimum gRNA length` (defined in the config file). As well as finding all gRNAs, these programs also produce many false positives, duplicates and overlapping alignments. In addition, no checks are made for an anchor, nor whether an alignment contains insertions or deletions in the mRNA. These checks are all done at a later step. These programs output to two files `maxi_alignments.txt` and `mini_alignments.txt`. 

#### Step 6

Find high quality canonical gRNAs with `hq_gRNAs()`.

To help in identifying inverted repeats and canonical gRNAs we first find the high quality canonical gRNAs; those that we can be quite sure are actual gRNAs. 

High quality gRNAs should be long (>= 40 nt) with a minimal number of mis-matches (<= 1), have a long anchor (>= 8 nt) of only Watson-Crick basepairs. These parameters are defined in the config file under `high quality gRNAs filter`.

On running the function `hq_gRNAs()`, the alignments from Step 5 will be loaded and the high quality gRNAs extracted. The number of high quality gRNAs will be printed and the number of minicircles that contain these gRNAs. The gRNA lengths are plotted as a histogram. The histogram should be right-skewed with lengths up into the 50s. Also plotted as a histogram are the 5' positions of the high quality gRNAs on the minicircles. These should be clustered into one or more peaks which represent the cassettes. 

If no, or very few, high quality gRNAs are found then the filtering parameters in the config file should be weakened and  `hq_gRNAs()` re-run.

The function outputs a fasta file of the HQ gRNAs plus flanking sequences upstream and downstream of the gRNA. These sequences are used in Step 7 to identify common motifs, such as the inverted repeats and gRNA initiation sequences. The number of nucleotides upstream and downstream from the 5' end of each HQ gRNA are given by the parameters `upstream` and `downstream` in the config file. If no inverted repeat motifs are found in Step 7 these values will need to be increased.

#### Step 7

Search for inverted repeats and the initiation sequence using `meme`. 

If Meme is installed then it can be run from `pipeline.py` directly. This will search for the common motifs and save the results in the directory `Work_files/Meme`. 

If Meme is not installed then use the https://meme-suite.org/meme/ server. Change the meme-suite parameters to search for at least 3 motifs with a minimum width of 5. This should pick up the initiation sequence just upstream of the HQ gRNAs. Once Meme has finished, download the file `meme.txt` and save it in the directory `Work_files/Meme`.


#### Step 8

Extract motifs from Meme output with the function `extract_motifs()`.

With a web-browser open the page `meme.html` either in `Work_files/Meme` or on the meme-suite server. 

Meme will likely find motifs that are too long. You need to tell the pipeline how long the motifs should be with the parameters `repeat length` and `initiation sequence length` in the config file. You also need to tell the pipeline how many nucleotides to trim off the lefthandside of each of Meme's motifs with the parameters `forward repeat left trim`, `initiation sequence left trim` and `reverse repeat left trim`.

Once you've set these parameters run `extract_motifs()`. This will produce histograms of the positions of the forward repeat, the initiation sequence and the reverse repeat in the HQ gRNA sequences saved in Step 6 (top row of histograms). It will also plot histograms of the distances from the forward repeat to the initiation sequence, the initiation sequence to the HQ gRNA and the forward repeat to the HQ gRNA. These should be used to make sure the positions of the motifs look okay. The histogram plot can be turned off by setting `plot feature histograms` to `no`.

The function `extract_motifs()` outputs two files: `motifs text file` which contains the regular expressions of the motifs and the nucleotide frequencies for each position in each motif; and `motifs pickle file` which is used in Step 9.

#### Step 9

Run `mO_scoring()` to create scoring vectors for inverted repeats and gRNAs (canonical and non-canonical) for each minicircle by using the nucleotide frequencies (also called biases) of these features from HQ gRNAs. 

The function needs to know the `expected number of cassettes` per minicircle. This is set in the config file. If not known, the expected number of cassettes can be guessed from the histogram of 5' positions of high quality gRNAs produced by `hq_gRNAs()`.

This function outputs `scores pickle file` for cassette identification in Step 10.

#### Step 10

Identify cassettes and label them with `identify_cassettes()`.

Cassettes are automatically identified using the scoring vectors from Step 9. The difficult part is assigning cassette labels to each cassette because the positions of cassettes can be quite variable between minicircles.

The function plots a histogram of the 5' ends of the forward repeats. The aim is to assign each cassette a label (e.g., I, II, III, etc) based on its position on it's minicircle. The config file has a dictionary parameter `cassette labels and limits` that requires the user to define a label and a closed set of positions on the minicircle. For example,

    cassette labels and limits: 
        I:   [100, 259]
        II:  [260, 409]
        III: [410, 469]
        IV:  [470, 703]
        V:   [653, 2000]

This says that cassettes whose 5' ends start between positions 100 and 259 inclusive should have label I. The ranges for some cassette labels are obvious from the histogram. However, the histograms of other cassettes may overlap making it difficult to assign a range to a label. The function assigns labels sequentially from 5' to 3' based on this dictionary. **Ranges may overlap** to account for the overlap in the histograms. For example, cassette labels IV and V overlap from positions 653 to 703.

If the function cannot assign a label to a cassette based on this dictionary it will output a warning and list all the cassettes on any offending minicircles. It is up to the user to change the dictionary and re-run the function until all cassettes can be assigned a label.

Sometimes it might be impossible to resolve all conflicts. This is usually because a cassette has been identified which is actually a false positive. In these cases the scoring vectors can be plotted for offending minicircles to help identify these false positive cassettes. List any such minicircles in the parameter `minicircles to plot`. For example, to plot the scoring vectors for "mO_097" do `minicircles to plot: [mO_097]` and re-run the function. 

The scoring plot shows the scores for the forward repeat (black) and the reverse repeat (red). Cassettes occur where these scores exceed the thresholds given as horizontal dashed lines of the same colour. A cassette is predicted to occur when a pair of black and red vertical lines appear together separated by roughly 100nt. The gRNA score is shown as a smooth blue line and gRNAs are predicted where peaks occur above the gRNA threshold shown as a dotted blue horizontal line. Above the scores are horizontal orange lines indicating the predicted positions of the cassettes, and blue crosses indicate the 5' ends of predicted gRNAs. A false positive cassette will have a pair of vertical black and red lines but a gRNA peak below the threshold. False positive cassettes can be manually removed by including their minicircle name and position in the the config file parameter `cassettes to drop` like so:

    cassettes to drop:
        - [mO_097, 120]

#### Step 11

Identify all canonical gRNAs by running `identify_all_gRNAs()`.

With cassettes identified in Step 10, it is now possible to identify all canonical gRNAs in cassettes and orphan gRNAs outside of cassettes. The filter parameters for gRNAs are given by the parameter `all gRNAs filter`.

Sometimes obvious false positive gRNAs come to light later that haven't been filtered out in this Step. These can be manually removed by listing their minicircle name and cassette label in parameter `false_positives` and the function re-run.

If transcriptomics is available for identifying expression status and initiation sequences then `have transcriptomics` should be set to "yes". Go to Step 12.

If transcriptomics is not available then the position of the putative initiation site, relative to the 3' end of the forward repeat, is the modal distance (calculated from the HQ gRNAs) between the forward repeat and the initiation sequence motifs as found by Meme in Steps 7 and 8. Anchors are identified and gRNAs are assigned to gRNA families. Human readable text files of all cassettes and canonical gRNAs are output to the `annotation directory` in files `cassettes text file` and `gRNAs text file` respectively. Go to Step 17.

#### Step 12

Assign small RNA transcripts to cassettes with `transcripts()`.

Using the sam file of transcripts aligned to minicircle sequences, these transcripts are now assigned to a cassette. 

A histogram of the position of the 5' ends of the templated sequences of the transcripts relative to the 3' end of the forward repeat is plotted to show the range of positions of the initiation site. The config file parameter `initiation site range` should be manually changed to this range for later use. 

#### Step 13

Find the probability, _p_, of a randomly selected transcript being in the initiation site under the null hypothesis that the cassette is not expressed by running `find_transcript_p()`. 

This function requires `initiation site range` from Step 12. The distribution of the position of the 5' ends of the templated sequences of the transcripts relative to the 3' end of the forward repeat is plotted on a log scale. The fit of a normal distribution curve is shown as a blue line. The transcripts within the initiation site are not included in this distribution. This data for this histogram is saved for later use in the config file parameter `transcript position distribution pickle file`.

The aim is to get a good fit of a normal distribution curve to the empirical distribution in order to calculate _p_.  The parameter in the config file `transcript position fit range` should be adjusted to trim the tails of the empirical distribution in order to obtain a good fit. 

The function outputs the polynomial coefficients of the fit, the estimated mean and standard of the normal curve and the probability _p_. The config file parameter `p` should be assigned the value of _p_.

#### Step 14

Find the statistic with the minimum variance to predict the end position of transcribed gRNA genes by running `predict_transcript_end_pos()`.

The transcripts are used to predict the transcribed end position of gRNA genes. As a cassette's aligned transcripts end at different positions on the minicircle we need a statistic that best represents the transcribed end of the cassette's gRNA gene.  This function outputs histograms of various statistics of the distance between the 3' end of gRNAs found by alignment and the 3' ends of the aligned transcripts. The statistic with the minimum variance is recommended by the function. The config file parameter `end position percentile` should be set to this value.

#### Step 15

Predict the expression status of sense and anti-sense strands of all cassettes by running `predict_expression()`.

Based on the value of _p_ and the number of sense and anti-sense transcripts aligning to each cassette and its initiation site we calculate the probability of each cassette's strand being expressed or not. 

The function outputs human readable expression information to `expression text file` and a pickle file for Step 16. 

#### Step 16

Add expression information to cassettes and canonical gRNAs and identify transcribed gRNA genes in each cassette by running `add_expression()`.

The expression information found in Step 15 is added to the information on cassettes and canonical gRNAs. If the `trim_to_init` value of the config file parameter `all gRNAs filter` is set to "yes" then any canonical gRNAs that begin upstream of the initiation sequence will have their 5' ends trimmed to start of the initiation sequence. In addition each canonical gRNA's anchor and gRNA family are found. 

Canonical and non-canonical gRNA genes in expressed cassettes are identified using the predicted end positions found in Step 14. 

Human readable text files of all cassettes, canonical gRNAs and genes are output to the `annotation directory` in files `cassettes text file`, `gRNAs text file` and `expressed genes text file` respectively. 

The pickle file `features_with_expression.pickle.gz` us output with all features and including minicircle and mRNA sequences for annotation output and analysis.

#### Step 17

Run `annotate_minicircles()` to output a genbank file of all the minicircles and their annotated features in `genbank directory/genbank text file`and alignments of canonical gRNAs to edited mRNA in the `alignments directory`.  There is one file per mRNA.

#### Step 18

Analysis. TODO.
