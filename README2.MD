The known CSB-1, 2 and 3 sequences as regular expressions are given in the config file. 

If there are no matches for CSB-1 in a minicircle's sequence this could because the minicircle has not been assembled completely, the regular expression is wrong or that a minor variant of the CSB-1 sequence exists. The program will output possible new variants that it finds at the start the minicircle sequences and then terminate. These variants can be added to the CSB-1 regular expression in the config file.

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