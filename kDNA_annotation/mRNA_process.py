import pathlib

from .common import load_config, get_directories, SeqIO

def get_unprocessed_mRNAs(unedited_file, edited_file):
    """ 
        Process fasta files of edited and unedited mRNA sequences and output
        files for minicircle annotataion

        The non-edited 5' and 3' ends of the edited and unedited sequences must
        match otherwise the processing gets stuck    
    """
    unedited = SeqIO.to_dict(SeqIO.parse(unedited_file, 'fasta'))
    edited   = SeqIO.to_dict(SeqIO.parse(edited_file,   'fasta'))
    mRNAs    = {}

    # for each mRNA
    for mRNA_name, mRNA in sorted(edited.items()):
        print(mRNA_name)
        # convert U's to T's in edited mRNAs
        edited_seq = str(mRNA.seq).replace('u', 'T').replace('U', 'T')

        # get mRNA name without version number for processing unedited sequences
        u_mRNA_name = mRNA_name.split('_')[0]
        unedited_seq = str(unedited[u_mRNA_name].seq).replace('u', 'T').replace('U', 'T')

        # by comparing edited and unedited sequences we can work out where 
        # insertions and deletions occur
        e_pos, u_pos = 0, 0
        deletion_list = []
        insertion_list = []
        unedited_list = []
        while e_pos < len(edited_seq):
            be = edited_seq[e_pos]
            bu = unedited_seq[u_pos]
            if bu == be:
                # edited and unedited bases are the same
                deletion_list.append('-')
                insertion_list.append(be)
                unedited_list.append(be)
                e_pos += 1
                u_pos += 1
            elif bu != be:
                # edited and unedited bases are not the same
                if be == 'T':
                    # if edited base is T then this is an insertion
                    deletion_list.append('-')
                    insertion_list.append('t')
                    unedited_list.append('-')
                    e_pos += 1
                elif bu == 'T':
                    # else if unedited base is T then this is a deletion
                    if deletion_list[-1] == '-':
                        deletion_list[-1] = 1
                    else:
                        deletion_list[-1] += 1
                    unedited_list.append(bu)
                    u_pos += 1
                else:
                    print("non-edited 5' or 3' sequences do not match in edited and unedited sequences")
                    print('Unedited sequence')
                    print(unedited_seq)
                    print('Edited sequence')
                    print(edited_seq)
                    print(f'Positions not matching: {u_pos}, {e_pos}')
                    print(f'Bases not matching: {bu}, {be}')
                    exit()

        constructed_deletions = ''.join([str(i) for i in deletion_list])
        constructed_edited_seq = ''.join(insertion_list)
        constructed_unedited_seq = ''.join(unedited_list)
        mRNAs[mRNA_name] = (constructed_edited_seq, constructed_deletions, constructed_unedited_seq)
    return mRNAs

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    in_dir, work_dir, _, meme_dir = get_directories(config)

    unedited_in_file  = f"{in_dir}/{config['unedited mRNA fasta infile']}"
    edited_in_file    = f"{in_dir}/{config['edited mRNA fasta infile']}"
    unedited_out_file = f"{work_dir}/{config['unedited mRNA fasta file']}"
    edited_out_file   = f"{work_dir}/{config['edited mRNA fasta file']}"
    edited_u_out_file = f"{work_dir}/{config['edited mRNA with u fasta file']}"
    edited_t_out_file = f"{work_dir}/{config['edited mRNA with t fasta file']}"
    deletion_out_file = f"{work_dir}/{config['deletions mRNA text file']}"

    pathlib.Path(work_dir).mkdir(parents=True, exist_ok=True) 
    pathlib.Path(meme_dir).mkdir(parents=True, exist_ok=True) 


    ####################################### PROCESS mRNA #########################################
    # read in unedited and edited mRNA sequences
    mRNAs = get_unprocessed_mRNAs(unedited_in_file, edited_in_file)


    ####################################### SAVE #################################################
    #output edited mRNA with all U's and u's as T's
    with open(edited_out_file, 'w') as f:
        for mRNA_name, seqs in mRNAs.items():
            f.write('>{}\n'.format(mRNA_name))
            f.write('{}\n'.format(seqs[0].replace('t', 'T')))

    #output positions of deletions
    with open(deletion_out_file, 'w') as f:
        for mRNA_name, seqs in mRNAs.items():
            f.write('>{}\n'.format(mRNA_name))
            f.write('{}\n'.format(seqs[1]))

    #output edited mRNA with t for an insertion
    with open(edited_t_out_file, 'w') as f:
        for mRNA_name, seqs in mRNAs.items():
            f.write('>{}\n'.format(mRNA_name))
            f.write('{}\n'.format(seqs[0]))

    with open(edited_u_out_file, 'w') as f:
        for mRNA_name, seqs in mRNAs.items():
            f.write('>{}\n'.format(mRNA_name))
            f.write('{}\n'.format(seqs[0].replace('t', 'u').replace('T', 'U')))

    with open(unedited_out_file, 'w') as f:
        for mRNA_name, seqs in mRNAs.items():
            f.write('>{}\n'.format(mRNA_name))
            f.write('{}\n'.format(seqs[2].replace('T', 'U').replace('-', '')))

