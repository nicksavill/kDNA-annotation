from .common import *

def get_motif_pos(seq, regex):
    match = re.search(regex, seq)
    if match:
        return match.start(0)
    else:
        return pd.NA

def realign_on_CSB1(minicircle):
    # If no CSB1 or already aligned do nothing
    if minicircle['CSB1_pos'] is pd.NA or minicircle['CSB1_pos'] == 0:
        return minicircle['seq']
    return minicircle['seq'][minicircle['CSB1_pos']:] + minicircle['seq'][:minicircle['CSB1_pos']]
        

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    in_dir, work_dir = get_directories(config)[:2]

    minicircle_file = f"{in_dir}/{config['minicircle fasta infile']}"
    maxicircle_file = f"{in_dir}/{config['maxicircle fasta infile']}"
    minicircle_clean_file = f"{work_dir}/{config['minicircle clean fasta file']}"
    maxicircle_clean_file = f"{work_dir}/{config['maxicircle clean fasta file']}"
    clean_log_file = f"{work_dir}/{config['clean log text file']}"


    ########################################## PARAMETERS #########################################
    CSB_regexes = config['CSB regexes']
    remove_non_CSB1 = config['remove minicircles without CSB1']


    ####################################### LOAD SEQUENCES #########################################
    minicircles = get_minicircles(minicircle_file)
    maxicircle  = get_maxicircle(maxicircle_file)

    # convert minicircles to a DataFrame for ease of use
    names = list(minicircles.keys())
    seqs = [str(m.seq) for m in minicircles.values()]
    description = [m.description for m in minicircles.values()]
    minicircles = pd.DataFrame({'mO_name':names, 'description':description, 'seq':seqs})


    ##################################### CHECKS ##################################################
    # if changed is True output output a warning 
    changed = False

    # check for duplicates
    dups = minicircles.duplicated(subset='seq')
    if dups.any():
        print('## The following sequences are duplicates. Keeping only the first instance of each duplicate.')
        print(minicircles[dups]['description'].to_string(index=False), end='\n\n')
        minicircles = minicircles[~dups].reset_index(drop=True)
        changed = True

    # check for CSB1
    # get length of CSB1 motif
    CSB1_len = len(re.sub('\[[ACGT]+\]', 'N', CSB_regexes['CSB1']))
    CSB1_missing = ~minicircles['seq'].str.contains(CSB_regexes['CSB1'])
    if CSB1_missing.any():
        print('## The following minicircles do not contain CSB1')
        print(minicircles[CSB1_missing]['description'].to_string(index=False), end='\n\n')
        if remove_non_CSB1:
            print('Removing minicircles without CSB1\n')
            minicircles = minicircles[~CSB1_missing]
            changed = True
        else:
            print('Possible CSB1 variants found in these minicircles. If valid, add them to configuration file and rerun.')
            variants = minicircles[CSB1_missing]['seq'].str[:CSB1_len]
            print(variants.unique(), end='\n\n')
            print('Exiting')
            exit()

    # check sequences are aligned on CSB1
    not_CSB1_aligned = ~minicircles['seq'].str.match(CSB_regexes['CSB1'])
    if not_CSB1_aligned.any():
        print('## The following minicircles are not aligned on CSB1')
        print(minicircles[not_CSB1_aligned]['description'].to_string(index=False), end='\n\n')
        print('Realigning, although this may not work if sequence contains multiple CSB1 motifs.\n')
        minicircles['CSB1_pos'] = minicircles['seq'].apply(get_motif_pos, args=(CSB_regexes['CSB1'],))
        minicircles['seq'] = minicircles.apply(realign_on_CSB1, axis=1)
        changed = True

    # check CSB3 is between 80 and 120nt downstream from CSB1
    minicircles['CSB3_pos'] = minicircles['seq'].apply(get_motif_pos, args=(CSB_regexes['CSB3'],))
    CSB3_wrong_pos = (80 < minicircles['CSB3_pos']) & (minicircles['CSB3_pos'] > 120)
    if CSB3_wrong_pos.any():
        minicircles['CSB1_count'] = minicircles['seq'].str.count(CSB_regexes['CSB1'])
        print('## The following minicircles have CSB3 in the wrong position')
        print(minicircles[CSB3_wrong_pos][['description', 'CSB3_pos', 'CSB1_count']].to_string(index=False), end='\n\n')
        print('This may be because the sequence contains multiple CSB1 motifs and is aligned on the wrong one.')
        print('Manually edit the sequences to align on correct CSB1-CSB3 pair and rerun.')
        print('Exiting')
        exit()

    # check naming convention
    invalid_name = ~minicircles['mO_name'].str.match('mO_\d{3}')
    if (invalid_name).any():
        print('## Minicircles have invalid names. Renaming all minicircles in ascending numerical order.')
        print(f'See {clean_log_file} for details.\n')
        # print(minicircles[invalid_name]['description'].to_string(index=False), end='\n\n')
        # move index into a column, index must be in consecutive ascending order
        # reset index to get ascending numbers into the Index
        minicircles = minicircles.reset_index(drop=True)
        # reset again to get the index into column "index"
        minicircles = minicircles.reset_index()
        # make new minicircle name
        minicircles['mO_name'] = minicircles['index'].apply(lambda x: f'mO_{x+1:>03d}')
        changed = True


    if changed:
        print('############################ WARNING #################################')
        print(f'Re-writing minicircle fasta file as {minicircle_clean_file}')
        print('Changes have been made to the minicircle sequence fasta file.')
        print('If small RNA transcripts have been aligned to the minicircle sequences.')
        print('this will have to be re-run with the cleaned minicircle sequences.')
        print('############################ WARNING #################################')
    else:
        print('All minicircles fine.')

################################################ SAVE #################################################
    with open(minicircle_clean_file, 'w') as f:
        for _, mO in minicircles.iterrows():
            f.write(f'>{mO["mO_name"]}\n')
            f.write(f'{mO["seq"]}\n')

    with open(maxicircle_clean_file, 'w') as f:
        f.write('>Maxicircle\n')
        f.write(f'{str(list(maxicircle.values())[0].seq)}\n')

    dataframe_out(minicircles[['mO_name', 'description']], clean_log_file, index=False)