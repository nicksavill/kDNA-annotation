import re
import pandas as pd
import matplotlib.pyplot as plt

from .common import *

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir = get_directories(config)[1]
    meme_dir = get_directories(config)[3]

    meme_file = f"{meme_dir}/meme.txt"
    motifs_text_file = f"{work_dir}/{config['motifs text file']}"
    motifs_pickle_file = f"{work_dir}/{config['motifs pickle file']}"


    ########################################## PARAMETERS #########################################
    up = config['upstream']
    down = config['downstream']
    repeat_len = config['repeat length']
    init_seq_len = config['initiation sequence length']
    plot_features = config['plot feature histograms']

    feature_names = ['forward repeat', 'init sequence', 'reverse repeat']

    new_meme_regex = re.compile('MEME-\d+\W+width')
    position_regex = re.compile('Start')
    nt_freq_regex  = re.compile('letter.probability matrix')
    regex_regex    = re.compile('regular expression')

    lengths = {
        'forward repeat':repeat_len,
        'init sequence':init_seq_len,
        'reverse repeat':repeat_len,
    }
    left_trim = {
        'forward repeat':config['forward repeat left trim'],
        'init sequence':config['initiation sequence left trim'],
        'reverse repeat':config['reverse repeat left trim'],
    }


    ################################### EXTRACT MOTIF INFORMATION ######################
    motif_positions = []
    motif_nt_freqs = {}
    motif_regex = {}
    this_motif, motif = None, None

    """
        Parse motif information from meme.txt
        states
            0 - searching for new motif
            1 - searching for motif positions around each HQ gRNA
            2 - recording motif positions around each HQ gRNA
            3 - searching for nt frequencies of motif
            4 - recording nt frequencies of motif
            5 - search and record regex of motif
    """
    state = 0
    with open(meme_file) as f:
        for line in f:
            if state == 0:
                # search for the next motif
                match = new_meme_regex.search(line)
                if match:
                    motif = line.split()[1]
                    motif_positions.append({'seq':[], motif:[]})
                    motif_nt_freqs[motif] = {'A':[], 'C':[], 'G':[], 'T':[]}
                    state = 1

            elif state == 1:
                # search for start of motif positions
                match = position_regex.search(line)
                if match:
                    next(f)
                    this_motif = motif_positions[-1]
                    state = 2

            elif state == 2:
                # recording motif positions
                try:
                    x = line.split()
                    seq = x[0]
                    start = int(x[1])
                    this_motif['seq'].append(seq)
                    this_motif[motif].append(start)
                except IndexError:
                    state = 3
                except ValueError:
                    print('Make sure to search for motifs on the positive sense strand only')
                    exit()
                    
            elif state == 3:
                # searching for nt frequencies
                match = nt_freq_regex.search(line)
                if match:
                    state = 4

            elif state == 4:
                try:
                    x = [float(i) for i in line.rstrip().split()]
                    motif_nt_freqs[motif]['A'].append(x[0])
                    motif_nt_freqs[motif]['C'].append(x[1])
                    motif_nt_freqs[motif]['G'].append(x[2])
                    motif_nt_freqs[motif]['T'].append(x[3])
                except ValueError:
                    state = 5

            elif state == 5:
                # search and record motif regular expression
                match = regex_regex.search(line)
                if match:
                    next(f)
                    motif_regex[motif] = next(f).rstrip()
                    # finished with this motif, start on the next one
                    state = 0


    # convert motif_positions list into a dataframe with each motif as a column and each HQ gRNA as a row
    motif_positions = pd.concat([pd.DataFrame(i).set_index('seq') for i in motif_positions], join='inner', axis=1)

    # calculate the median positions of each motif in order to convert motif seq into a feature name
    medians = motif_positions.median().sort_values()
    # dictionary of motif seq (key) and its corresponding feature name
    convert = dict([(i, j) for i, j in zip(medians.index, feature_names)])

    # convert motif seqs to feature names in motif_positions dataframe and order columns
    motif_positions = motif_positions.rename(columns=convert)
    motif_positions = motif_positions[feature_names]

    # convert motif nt frequency and regex keys from motif seq to feature name
    for motif, feature in convert.items():
        motif_nt_freqs[feature] = motif_nt_freqs.pop(motif)
        motif_regex[feature] = motif_regex.pop(motif)

    # trim unwanted nt frequencies
    for motif, nt_freqs in motif_nt_freqs.items():
        for base in nt_freqs:
            s = left_trim[motif]
            e = s+lengths[motif]
            nt_freqs[base] = nt_freqs[base][s:e]

    # trim left end of regular expression
    for motif, regex in motif_regex.items():
        r = f'(\[[ACGT]+\]|[ACGT]){{{left_trim[motif]}}}((\[[ACGT]+\]|[ACGT]){{{lengths[motif]}}})'
        match = re.match(r, regex)
        if match:
            motif_regex[motif] = match.group(2)
        else:
            print('Something wrong with motif regular expression for {motif} in meme.txt')

    # check repeats sit within the sequences flanking the HQ gRNA
    # if not get user to change upstream and downstream and rerun analysis
    motif_positions['forward repeat'] += left_trim['forward repeat']
    reverse_end = motif_positions['reverse repeat'] + left_trim['reverse repeat'] + repeat_len
    x = reverse_end.max() -(up+down)
    if x >= 0:
        print("*** WARNING ***")
        print(f"The 3' end of the reverse repeat has been truncated by up to {x:.0f} nt on some high quality gRNAs.")
        print(f"Increase 'upstream'  by at least {x:.0f} nt in config.yaml, and rerun hq_gRNA.py and meme again.")

    if plot_features:
        # check that features are where we expect them to be
        motif_positions['forward end to gRNA'] = up - repeat_len - motif_positions['forward repeat']
        motif_positions['init to gRNA'] = up - motif_positions['init sequence']
        motif_positions['forward end to init'] = motif_positions['init sequence'] - motif_positions['forward repeat'] - repeat_len

        fig, axes = plt.subplots(2, 3, figsize=(3*6.4, 2*4.8))
        axes[0, 0].hist(motif_positions['forward repeat'])
        axes[0, 1].hist(motif_positions['init sequence'])
        axes[0, 2].hist(motif_positions['reverse repeat'])
        axes[1, 0].hist(motif_positions['forward end to init'])
        axes[1, 1].hist(motif_positions['init to gRNA'])
        axes[1, 2].hist(motif_positions['forward end to gRNA'])

        axes[0, 0].set_xlabel("Position of 5' end of forward repeat")
        axes[0, 1].set_xlabel("Position of 5' end of initiation sequence")
        axes[0, 2].set_xlabel("Position of 5' end of reverse repeat")
        axes[1, 0].set_xlabel("Distance from 3' end of forward repeat\nto 5' end of initiation sequence")
        axes[1, 1].set_xlabel("Distance from 5' end of initiation sequence\nto 5' end of HQ gRNA")
        axes[1, 2].set_xlabel("Distance from 3' end of forward repeat\nto 5' end of HQ gRNA")
        fig.suptitle('Position of features in sequences upstream and downstream of high quality gRNAs')
        plt.show()


    ##################################################### SAVE ###########################################################

    with open(motifs_text_file, 'w') as f:
        f.write('Motif regular expressions of HQ gRNAs\n')
        
        for motif, regex in motif_regex.items():
            f.write(f'{motif}:\t{regex}\n')
        f.write('\n')
        
        f.write('Motif nucleotide frequencies of HQ gRNAs\n')

        for motif, nt_freq in motif_nt_freqs.items():
            f.write(f'{motif}\n')
            f.write('       A      C      G      T\n')
            for i in range(lengths[motif]):
                f.write(f'{i:< 3d}  ')
                for b in nt_freq:
                    f.write(f'{nt_freq[b][i]:.3f}  ')
                f.write('\n')
            f.write('\n')

    pickle_save([motif_positions, motif_regex, motif_nt_freqs], motifs_pickle_file)
