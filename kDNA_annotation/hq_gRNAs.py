"""
Identify high quality gRNA found by alignment (this might include orphans)

Plot positions of high quality gRNAs identified by alignment to edited mRNA
to get an idea of their positions on minicircles

Output gRNA sequences including upstream and downstream flanking sequences to search for motifs using meme:
    1. forward and reverse inverted repeats
    2. initiation sequences

Output is saved as
    hq_gRNAs.pickle - pickle file for later import
    hq_gRNAs.txt    - text file of dataframe for manual checking
    hq_gRNAs.fasta  - fasta file for MEME analysis
"""

import matplotlib.pyplot as plt
import re
from Bio.Seq import reverse_complement
from collections import OrderedDict

from .common import *

def extract_alignment(alignment, gRNAs, mRNAs, circle_name, minicircles, filter):
    mRNA_name, mRNA_align_o, gRNA_align_o, strand = alignment.split()
    circle_name = circle_name.split()[0]
    if mRNA_name in mRNAs:
        mRNA_record = mRNAs[mRNA_name]
    else:
        return

    dd = {'c':'coding', 't':'template'}
    strand = dd[strand]

    # get the pairing for this alignment
    mRNA_seq = mRNA_align_o.replace('T', 'U')
    gRNA_seq = complement(gRNA_align_o)
    if gRNA_seq == None:
        return
    else:
        gRNA_seq = gRNA_seq.replace('T', 'U')
    pairing  = pairs(mRNA_seq, gRNA_seq)

    # find 5'-most WC anchor region
    min_anchor_length = filter['min_anchor_length']
    min_gRNA_length = filter['min_gRNA_length']
    anchor_seq_regex = re.compile(filter['anchor_seq'])

    anchors = anchor_seq_regex.finditer(pairing)
    try:
        # get the longest 5'-most anchor
        anchor = [a.group(1) for a in anchors if len(a.group(1)) >= min_anchor_length][-1]
        anchor_pos = pairing.rfind(anchor)
        length = anchor_pos+len(anchor)
    except IndexError:
        return
    if length < min_gRNA_length:
        return

    # trim alignment so that it starts at the anchor
    mRNA_align = mRNA_align_o[:length]
    gRNA_align = gRNA_align_o[:length]
    pairing = pairing[:length]

    # find position on mRNA of gRNA and check edits occur there
    mRNA_start = mRNA_record['DNA_seq'].index(mRNA_align)
    if 't' not in mRNA_record['edits'][mRNA_start:mRNA_start+length]:
        return

    # find start of gRNA on minicircle and assign cassette_pos or orphan or Maxi
    if strand == 'template':
        start = str(minicircles[circle_name].seq).index(gRNA_align)
    else:
        start = str(minicircles[circle_name].seq).index(reverse_complement(gRNA_align))

    # saved alignment information
    gRNA = OrderedDict()
    gRNA['mO_name']      = circle_name
    gRNA['strand']       = strand
    gRNA['length']       = length
    gRNA['circle_start'] = start
    gRNA['circle_end']   = start+length
    gRNA['mRNA_name']    = mRNA_name
    gRNA['product']      = mRNA_name.split('_')[0]
    gRNA['mRNA_start']   = mRNA_start
    gRNA['mRNA_end']     = mRNA_start+length
    gRNA['mRNA_seq']     = mRNA_record['edits'][mRNA_start:mRNA_start+length].replace('t', 'u').replace('T', 'U')
    gRNA['gRNA_seq']     = complement(gRNA_align).replace('T', 'U')
    gRNA['pairing']      = pairing
    gRNA['mismatches']   = gRNA['pairing'].count('.')
    gRNA['anchor_len']   = len(anchor_seq_regex.match(pairing[::-1]).group(1))
    gRNAs.append(gRNA)

def get_hq_gRNAs(mini_align_file, minicircles, mRNAs, filter):
    # read in the alignments from minicircle to edited mRNA
    # these files contain many duplicates, overlapping alignments, low quality false positives
    gRNAs = []
    with open(mini_align_file) as f:
        for name in f:
            extract_alignment(next(f), gRNAs, mRNAs, name[1:-1], minicircles, filter)

    # convert alignments to dataframe for filtering
    gRNAs = pd.DataFrame(gRNAs)
    if len(gRNAs) == 0:
        print('No high quality gRNAs found')
        exit()
    print(f'from alignment = {len(gRNAs)}')

    # filter out all duplicate alignments
    gRNAs = gRNAs.drop_duplicates()
    print(f'after drop duplicates = {len(gRNAs)}')

    # if >1 alignment begins or ends at same mapped position on an mRNA remove the smaller ones  
    # have to do this otherwise collapsing doesn't work if >1 gRNA in a cassette edits different parts of an mRNA
    gRNAs = gRNAs.groupby(['mO_name', 'strand', 'mRNA_name', 'circle_start', 'mRNA_end']).apply(drop_smaller).reset_index(drop=True)
    gRNAs = gRNAs.groupby(['mO_name', 'strand', 'mRNA_name', 'circle_start', 'mRNA_start']).apply(drop_smaller).reset_index(drop=True)
    gRNAs = gRNAs.groupby(['mO_name', 'strand', 'mRNA_name', 'circle_end', 'mRNA_start']).apply(drop_smaller).reset_index(drop=True)
    gRNAs = gRNAs.groupby(['mO_name', 'strand', 'mRNA_name', 'circle_end', 'mRNA_end']).apply(drop_smaller).reset_index(drop=True)
    print(f'after drop smaller sub-alignnments = {len(gRNAs)}')

    # trim the mismatches at the left-most end of the alignment
    gRNAs = gRNAs.apply(trim_mismatches, axis=1)

    # filter out low quality alignments, step 1
    minl = gRNAs['length'] < filter['min_gRNA_length']
    mina = gRNAs['anchor_len'] < filter['min_anchor_length']
    maxm = gRNAs['mismatches'] > filter['max_mismatches']
    gRNAs = gRNAs[~(minl | mina | maxm)]

    # # collapse overlapping alignments of same mRNA and same minicircle into a longer alignment
    gRNAs = gRNAs.groupby(['mO_name', 'strand', 'mRNA_name']).apply(collapse, *(minicircles, mRNAs), **{'post_cassette_label':False}).reset_index(drop=True)
    print(f'after collapsing = {len(gRNAs)}')

    # # trim a second time as some collapsed alignments extend too far
    gRNAs = gRNAs.apply(trim_mismatches, axis=1)
    
    # filter out all duplicate alignments again (this prevents double counting gRNAs of multiple mRNA versions)
    gRNAs = gRNAs.drop_duplicates(subset=['mO_name', 'strand', 'pairing'])
    print(f'after drop duplicates = {len(gRNAs)}\n')

    print('Number of HQ gRNAs on each strand')
    print(gRNAs['strand'].value_counts())
    major_strand = gRNAs['strand'].mode().loc[0]
    gRNAs = gRNAs[gRNAs['strand'] == major_strand]
    print(f'Only using gRNAs on the {major_strand} strand\n')

    print(f'Total high quality gRNAs found = {len(gRNAs)}')
    return gRNAs

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir = get_directories(config)[1]

    minicircle_file      = f"{work_dir}/{config['minicircle clean fasta file']}"
    mini_align_file      = f"{work_dir}/{config['minicircle alignments file']}"
    edited_mRNA_t_file   = f"{work_dir}/{config['edited mRNA with t fasta file']}"
    deletion_mRNA_file   = f"{work_dir}/{config['deletions mRNA text file']}"
    hq_gRNAs_pickle_file = f"{work_dir}/{config['high quality gRNAs pickle file']}"
    hq_gRNAs_text_file   = f"{work_dir}/{config['high quality gRNAs text file']}"
    hq_gRNAs_fasta_file  = f"{work_dir}/{config['high quality gRNAs fasta file']}"


    ########################################## PARAMETERS #########################################
    # parameters defining high quality gRNAs
    filter = config['high quality gRNAs filter']

    # number of nt upstream and downstream of start of gRNA to output for motif calling and nt bias scoring 
    up = config['upstream']  
    down = config['downstream']  


    ####################################### LOAD SEQUENCES #########################################
    minicircles = get_minicircles(minicircle_file)
    mRNAs = get_mRNAs(edited_mRNA_t_file, deletion_mRNA_file)
    hq_gRNAs = get_hq_gRNAs(mini_align_file, minicircles, mRNAs, filter)

    # get sequence around gRNA for motif searching
    hq_gRNAs['sequence'] = hq_gRNAs.apply(get_sequence, args=(up, down, minicircles), axis=1)

    # proportion of minicircles on which a high quality gRNA has been identified
    print(f"{len(hq_gRNAs['mO_name'].unique())} of {len(minicircles)} minicircles represented")


    ####################################### SAVE HQ gRNAS #########################################
    # save for later use for nucleotide bias scoring
    pickle_save(hq_gRNAs, hq_gRNAs_pickle_file)

    # output hq_gRNAs to file to check
    dataframe_out(hq_gRNAs, hq_gRNAs_text_file)

    # output sequences of genes and their flanks for motif identification with meme
    # OR use http://meme-suite.org/

    with open(hq_gRNAs_fasta_file, 'w') as f:
        for _, gRNA in hq_gRNAs.iterrows():
            f.write(f'>{gRNA["mO_name"]}_{gRNA["circle_start"]}_{gRNA["mRNA_name"]}\n')
            f.write(f'{gRNA["sequence"]}\n')


    ################################################ CHECK #########################################
    # check gRNA lengths are okay and their positions on minicircles
    _, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].hist(hq_gRNAs['length'])
    axes[0].set_xlabel('gRNA length')
    axes[1].hist(hq_gRNAs['circle_start'], bins=range(0, hq_gRNAs['circle_start'].max()+1, 10))
    axes[1].set_xlabel('gRNA start positions in minicircles')
    plt.show()

