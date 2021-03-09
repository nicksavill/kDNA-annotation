from collections import OrderedDict
import re
import gzip
import matplotlib.pyplot as plt

from .common import *

def extract_alignments(cassettes, orphans, transcripts_alignment_file):
    """ extract all transcripts from SAM file and associate with cassettes and orphans
        Split transcript into templated and non-templated regions and
        find relative position in cassette to end of forward repeat. 
        Save pandas dataframe of all alignments
    """
    cigar_regex = re.compile('([0-9]+)([MIDNSHPX=])')
    strand_code = {'0':'coding', '16':'template'}

    def get_alignment(line):
        x = line.rstrip().split('\t')
        return [strand_code[x[1]]]+[x[2]]+[int(x[3])-1]+[int(x[4])]+[x[5]]+[x[9]]

    def edit_seq(seq, cigar, strand):
        if 'I' in cigar or 'D' in cigar:
            return None, None, None
        cigar_list = list(cigar_regex.findall(cigar))
        n5, n3 = 0, len(seq)

        # n5 is the length of the 5' non-coding seq (and the start of the coding seq)
        # n3 is the length of the 5' non-coding seq plus the length of coding seq

        # template strand transcripts in the bam file are 3' to 5' (left to right) 
        # but the complement of their actual sequence. We keep the same orientation (3' to 5')
        # but take the complement so they have the sequence they naturally have
        if strand == 'template':
            seq = complement(seq)

        if cigar_list[0][1] == 'S':
            n5 = int(cigar_list[0][0])
            seq = seq[:n5].lower()+seq[n5:]

        if cigar_list[-1][1] == 'S':
            n3 -= int(cigar_list[-1][0])
            seq = seq[:n3]+seq[n3:].lower()

        return seq, n5, n3

    columns = (
        'mO_name',
        'cassette_label',
        'strand',
        'mO_pos',
        'rel_pos',
        'end_pos',
        'non_coding_5p_seq',
        'coding_seq',
        'non_coding_3p_seq'
    )

    # create a dictionary of cassettes for fast lookup
    cas_dict = {}
    for _, c in cassettes.iterrows():
        mO = c['mO_name']
        if mO not in cas_dict:
            cas_dict[mO] = []
        cas_dict[mO].append((c['forward_end'], c['reverse_start'], c['cassette_label']))

    # allow a 20 nt window either side of orphan genes identified by alignment for associated transcripts
    window = 20
    orphan_dict = {}
    for _, c in orphans.iterrows():
        mO = c['mO_name']
        if mO not in orphan_dict:
            orphan_dict[mO] = c['circle_start']-window, c['circle_end']+window, c['cassette_label']

    df_data = OrderedDict([(c, []) for c in columns])

    """
        Currently this reads a gizpped SAM file without @SQ headers.
        This needs re-writing to use pysam to be able to read in transcripts 
        located in specific cassettes from indexed BAM files.
    """
    with gzip.open(transcripts_alignment_file, 'rt') as o:
        for line in o:
            strand, mO_name, pos, _, cigar, seq = get_alignment(line)
            if mO_name in cas_dict:
                if mO_name.startswith('mO') and pos > 0:
                    # first check if transcript is in a cassette
                    for c in cas_dict[mO_name]:
                        if c[0] < pos and pos+len(seq) < c[1]:
                            break
                    else:
                        # if not check that it aligns with an orphan gRNA
                        if mO_name not in orphan_dict:
                            continue
                        c = orphan_dict[mO_name]
                        if pos < c[0] or c[1] < pos+len(seq):
                            continue

                    seq, n5, n3 = edit_seq(seq, cigar, strand)

                    if seq is not None:
                        l = n3-n5
                        df_data['mO_name'].append(mO_name)
                        df_data['cassette_label'].append(c[2])
                        df_data['strand'].append(strand)
                        df_data['mO_pos'].append(pos)
                        df_data['non_coding_5p_seq'].append(seq[:n5])
                        df_data['coding_seq'].append(seq[n5:n3])
                        df_data['non_coding_3p_seq'].append(seq[n3:])

                        if strand == 'coding':
                            if c[2] == 'Orphan':
                                relpos = pos-(c[0]+window)
                            else:
                                relpos = pos-c[0]
                        else:
                            if c[2] == 'Orphan':
                                relpos = c[1]-window-(pos+l)
                            else:
                                relpos = c[1]-(pos+l)
                        df_data['rel_pos'].append(relpos)
                        # the position of the end of the transcript relative to 3' end of forward repeat
                        # or for orphan the start of the aligned gRNA
                        df_data['end_pos'].append(relpos+l)
    return pd.DataFrame(df_data)

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    in_dir, work_dir = get_directories(config)[:2]

    transcripts_alignment_file = f"{in_dir}/{config['transcriptomics infile']}"
    transcripts_file           = f"{work_dir}/{config['transcripts pickle file']}"
    features_file              = f"{work_dir}/{config['features pickle file']}"


    ########################################## PARAMETERS #########################################

    ############################################### LOAD ###########################################
    cassettes, gRNAs = gzip_pickle_load(features_file)[-2:]


    ################################# EXTRACT TRANSCRIPTS #########################################
    orphans = gRNAs.query('cassette_label == "Orphan"')

    # extract transcripts from sam file and associate to cassettes and orphans
    transcripts = extract_alignments(cassettes, orphans, transcripts_alignment_file)


    ########################################### SAVE #############################################
    gzip_pickle_save(transcripts, transcripts_file)


    ######################### PLOT TO FIND POSITION AND RANGE OF INITIATION SITE #################
    plt.hist(transcripts['rel_pos'], bins=range(transcripts['rel_pos'].max()))
    plt.xlabel('Transcript position relative to forward repeat')
    plt.ylabel('Frequency')
    plt.title("Distribution of transcript position relative to 3' end of forward repeat")
    plt.show()


    # TODO ADD ALIGNMENTS HERE


