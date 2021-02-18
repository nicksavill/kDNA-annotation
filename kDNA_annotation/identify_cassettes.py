"""
User-defined closed intervals on the minicircle for each cassette label. These intervals can overlap.
Intervals are user-defined based on the histogram of forward repeat start positions.
These may not work initially because of invalid cassettes or cassettes not quite sitting in the correct
positions on minicircles. If the algorithm cannot assign cassette labels correctly it outputs which are 
wrong and all the cassettes in the offending minicircle. Based on this, the intervals can be adjusted.  


Output is saved in
    cassettes.pickle
    cassettes.txt
"""

import matplotlib.pyplot as plt
import numpy as np

from .common import *

def get_cassettes(
        mO_forward_scores,
        mO_reverse_scores,
        gRNA_peaks,
        min_forward_score,
        _,
        min_reverse_score,
        forward_search_min,
        forward_search_max,
        reverse_search_min,
        reverse_search_max,
        repeat_len,
        minicircles
    ):    
    cassettes = []
    for _, group in gRNA_peaks.groupby('mO_name'):
        mO_name = group.iloc[0, 0]

        cassettes_pos = []
        q_peaks = group[group['quality'] == True]
        forward_score = mO_forward_scores[mO_name]
        reverse_score = mO_reverse_scores[mO_name]
        for _, peak in q_peaks.iterrows():
            x = peak['position']
            # calculate search region for forward repeat and find position of max score
            x1 = max(0, x + forward_search_min)
            x2 = x + forward_search_max
            fpos = np.argmax(forward_score[x1:x2]) + x1
            # calculate search region for reverse repeat and find position of max score
            x1 = fpos + reverse_search_min
            x2 = fpos + reverse_search_max
            rpos = np.argmax(reverse_score[x1:x2]) + x1

            # accept as cassette if either repeat score is above min repeat score 
            # and this cassette des not overlap existing one
            if forward_score[fpos] > min_forward_score or reverse_score[rpos] > min_reverse_score:
                for c in cassettes_pos:
                    if fpos < c[1] and rpos > c[0]:
                        break
                else:
                    cassettes_pos.append((fpos, rpos))

        # find high scoring repeats not within cassettes that have not yet been associated with a cassette 
        # find positions of probable repeats
        f_indexes = np.where(forward_score > min_forward_score)[0]
        r_indexes = np.where(reverse_score > min_reverse_score)[0]
        # and remove any that have already been assigned to a cassette
        for c in cassettes_pos:
            f_indexes = f_indexes[~((f_indexes >= c[0]) & (f_indexes <= c[1]))]
            r_indexes = r_indexes[~((r_indexes >= c[0]) & (r_indexes <= c[1]))]
        # reverse scores at remaining un-associated reverse repeats

        for fpos in f_indexes:
            # search window for reverse repeat from the position of the forward repeat
            x1 = fpos + reverse_search_min
            x2 = fpos + reverse_search_max
            # unassociated reverse repeats within range of forward repeat
            r_indexes_window = r_indexes[(r_indexes >= x1) & (r_indexes <= x2)]

            if len(r_indexes_window) > 0:
                idx = np.argmax(reverse_score[r_indexes_window])
                rpos = r_indexes_window[idx]
                # ignore if overlapping with an existing cassette
                for c in cassettes_pos:
                    if fpos < c[1] and rpos > c[0]:
                        break
                else:
                    cassettes_pos.append((fpos, rpos))

        for c in cassettes_pos:
            new_c = {}
            new_c['mO_name'] = mO_name
            new_c['forward_start'] = c[0]
            new_c['forward_end'] = c[0]+repeat_len
            new_c['forward_seq'] = str(minicircles[mO_name].seq)[c[0]:c[0]+repeat_len]
            new_c['reverse_start'] = c[1]
            new_c['reverse_end'] = c[1]+repeat_len
            new_c['reverse_seq'] = str(minicircles[mO_name].seq)[c[1]:c[1]+repeat_len]
            cassettes.append(new_c)

    cassettes = pd.DataFrame(cassettes)
    return cassettes.sort_values(['mO_name', 'forward_start']).reset_index(drop=True)

def get_cassettes_old(
        mO_forward_scores,
        mO_reverse_scores,
        gRNA_peaks,
        min_forward_score,
        _,
        min_reverse_score,
        up,
        forward_search_min,
        forward_search_max,
        reverse_search_min,
        reverse_search_max,
        repeat_len,
        minicircles
    ):    
    cassettes = []
    for _, group in gRNA_peaks.groupby('mO_name'):
        mO_name = group.iloc[0, 0]

        cassettes_pos = []
        q_peaks = group[group['quality'] == True]
        forward_score = mO_forward_scores[mO_name]
        reverse_score = mO_reverse_scores[mO_name]
        for _, peak in q_peaks.iterrows():
            x = peak['position']
            # calculate search region for forward repeat and find position of max score
            x1 = max(0, x - up + forward_search_min)
            x2 = x - up + forward_search_max
            fpos = np.argmax(forward_score[x1:x2]) + x1
            # calculate search region for reverse repeat and find position of max score
            x1 = max(0, x - up + reverse_search_min)
            x2 = x - up + reverse_search_max
            rpos = np.argmax(reverse_score[x1:x2]) + x1

            # accept as cassette if either repeat score is above min repeat score 
            # and this cassette des not overlap existing one
            if forward_score[fpos] > min_forward_score or reverse_score[rpos] > min_reverse_score:
                for c in cassettes_pos:
                    if fpos < c[1] and rpos > c[0]:
                        break
                else:
                    cassettes_pos.append((fpos, rpos))

        # find high scoring repeats not within cassettes that have not yet been associated with a cassette 
        # find positions of probable repeats
        f_indexes = np.where(forward_score > min_forward_score)[0]
        r_indexes = np.where(reverse_score > min_reverse_score)[0]
        # and remove any that have already been assigned to a cassette
        for c in cassettes_pos:
            f_indexes = f_indexes[~((f_indexes >= c[0]) & (f_indexes <= c[1]))]
            r_indexes = r_indexes[~((r_indexes >= c[0]) & (r_indexes <= c[1]))]
        # reverse scores at remaining un-associated reverse repeats

        for fpos in f_indexes:
            # search window for reverse repeat from the position of the forward repeat
            x1 = fpos+reverse_search_min-forward_search_max
            x2 = fpos+reverse_search_max-forward_search_min
            # unassociated reverse repeats within range of forward repeat
            r_indexes_window = r_indexes[(r_indexes >= x1) & (r_indexes <= x2)]

            if len(r_indexes_window) > 0:
                idx = np.argmax(reverse_score[r_indexes_window])
                rpos = r_indexes_window[idx]
                # ignore if overlapping with an existing cassette
                for c in cassettes_pos:
                    if fpos < c[1] and rpos > c[0]:
                        break
                else:
                    cassettes_pos.append((fpos, rpos))

        for c in cassettes_pos:
            new_c = {}
            new_c['mO_name'] = mO_name
            new_c['forward_start'] = c[0]
            new_c['forward_end'] = c[0]+repeat_len
            new_c['forward_seq'] = str(minicircles[mO_name].seq)[c[0]:c[0]+repeat_len]
            new_c['reverse_start'] = c[1]
            new_c['reverse_end'] = c[1]+repeat_len
            new_c['reverse_seq'] = str(minicircles[mO_name].seq)[c[1]:c[1]+repeat_len]
            cassettes.append(new_c)

    cassettes = pd.DataFrame(cassettes)
    return cassettes.sort_values(['mO_name', 'forward_start']).reset_index(drop=True)

def drop_invalid_cassettes(cassettes, cassettes_to_drop):
    """ remove cassettes in cassettes_to_drop """
    y = cassettes[['mO_name', 'forward_start']].apply(lambda x: list(x) in cassettes_to_drop, axis=1)
    return cassettes[~y].reset_index(drop=True)

def get_cassette_label(mO_cassettes, cas_labels):
    """ Try to assign cassette labels in cas_labels to the cassettes in the minicircle
        based on the start positions of the forward repeats """

    # An iterator of cassette labels as defined in cas_labels
    label_values = iter(cas_labels.keys())
    # The list of positions assigned to each cassette
    cassettes_label_list = []

    # loop through cassettes and try to assign a label based on cas_labels
    for forward_start in mO_cassettes:
        # get the first cassette label
        label = next(label_values)
        # search through label positions to find the first that this cassette sits in
        while forward_start > cas_labels[label][1]:
            label = next(label_values)
        # check that the cassette is in the label interval
        if forward_start >= cas_labels[label][0]:
            # if it is assign the label to this cassette
            cassettes_label_list.append(label)
        else:
            # otherwise the cassette is not in the label interval. Output a warning and return None
            print(f"cassette {mO_cassettes.name} {forward_start} should not be in position {label}")
            return None
    return cassettes_label_list

def plot(
        mO,
        minicircles,
        cassettes,
        mO_gRNA_scores,
        mO_forward_scores,
        mO_reverse_scores,
        gRNA_peaks,
        min_gRNA_score,
        min_forward_score,
        min_reverse_score,
        up
    ):
    if isinstance(mO, int):
        mO_name = f'mO_{mO:>03d}'
    else:
        mO_name = mO

    if mO_name not in minicircles:
        print(f'{mO_name} is not a valid minicircle name')
        exit()

    scores = mO_gRNA_scores.query('mO_name == @mO_name')['smoothed'].array[0]
    peaks = gRNA_peaks.query('mO_name == @mO_name and quality==True')['position'].array

    forward_score = mO_forward_scores[mO_name]
    reverse_score = mO_reverse_scores[mO_name]
    f_indexes = np.where(forward_score > min_forward_score)[0]
    r_indexes = np.where(reverse_score > min_reverse_score)[0]

    _, axes = plt.subplots(1, 1, figsize=(14, 10))
    axes.plot(forward_score, color='k', label="5' end forward repeat score")
    axes.plot(reverse_score, color='r', label="5' end reverse repeat score")
    axes.plot(np.arange(len(scores))+up, scores, color='b', label="5' end gRNA score")
    for _, c in cassettes.query('mO_name == @mO_name').iterrows():
        axes.plot([c['forward_start'], c['reverse_start']], [0, 0], color='Orange')
    axes.scatter(peaks, np.zeros_like(peaks), color='b', marker='x', label="peak 5' end gRNA score")
    axes.scatter(f_indexes, np.zeros_like(f_indexes), color='k')
    axes.scatter(r_indexes, np.zeros_like(r_indexes), color='r')
    axes.axhline(min_forward_score, color='k', ls='--', label=f'forward repeat threshold {min_forward_score:.1f}')
    axes.axhline(min_reverse_score, color='r', ls='--', label=f'reverse repeat threshold {min_reverse_score:.1f}')
    axes.axhline(min_gRNA_score, color='b', ls=':', label=f'gRNA threshold {min_gRNA_score:.1f}')
    axes.set_title(mO_name)
    # axes.legend()
    plt.show()

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir = get_directories(config)[1]

    minicircle_file       = f"{work_dir}/{config['minicircle clean fasta file']}"
    scores_file           = f"{work_dir}/{config['scores pickle file']}"
    cassettes_pickle_file = f"{work_dir}/{config['cassettes pickle file']}"
    cassettes_text_file   = f"{work_dir}/{config['cassettes text file']}"
    motifs_pickle_file    = f"{work_dir}/{config['motifs pickle file']}"


    ########################################## PARAMETERS #########################################
    # whether to plot the histogram of cassette start position
    plot_hist = config['plot cassette histogram']

    # which minicircles to plot
    plot_mOs = config['minicircles to plot']

    # the suggested cassette labels and their limits
    cas_labels = config['cassette labels and limits']

    # invalid cassettes to drop
    cassettes_to_drop = config['cassettes to drop']

    # some wiggle room in the search region for repeats
    wiggle = config['wiggle']

    repeat_len = config['repeat length']
    up = config['upstream']


    ####################################### LOAD SEQUENCES AND SCORES #############################
    minicircles = get_minicircles(minicircle_file)
    mO_scores = pickle_load(scores_file)
    motif_positions = pickle_load(motifs_pickle_file)[0]

    max_len = max([len(m.seq) for m in minicircles.values()])


    ############## CALCULATE SEARCH REGION FOR FORWARD AND REVERSE REPEATS #######################
    cassette_length = motif_positions['reverse repeat'] - motif_positions['forward repeat']
    forward_search_min = motif_positions['forward repeat'].min() - up - wiggle
    forward_search_max = motif_positions['forward repeat'].max() - up + wiggle
    reverse_search_min = cassette_length.min() - wiggle
    reverse_search_max = cassette_length.max() + wiggle + 1
    print(f"forward search region relative to 5' end of gRNA: {forward_search_min} - {forward_search_max}")
    print(f"reverse search region relative to 5' end of forward repeat: {reverse_search_min} - {reverse_search_max}")


    ################################### IDENTIFY CASSETTES #######################################
    info = [forward_search_min, forward_search_max, reverse_search_min, reverse_search_max, repeat_len]
    cassettes = get_cassettes(*(mO_scores[1:]+info), minicircles)
    print(f'{len(cassettes)} found')


    ##################### PLOT CASSETTE POSITIONS ON MINICIRCLES TO AID LABELLING #################
    if plot_hist:
        plt.hist(cassettes['forward_start'], bins=range(max_len))
        plt.show()


    ################################### PLOT MINICIRCLE SCORES ######################################
    for mO in plot_mOs:
        plot(mO, minicircles, cassettes, *(mO_scores+[up]))


    ############################## DROP USED-DEFINED INVALID CASSETTES ##################################
    cassettes = drop_invalid_cassettes(cassettes, cassettes_to_drop)


    ################## ASSIGN CASSETTE POSITION ACCORDING TO THE USER-DEFINED CAS_LABELS ###############
    # get each cassette's position on the minicircle given by the used-defined dictionary cas_labels
    print(f'In config.yaml "cassette labels and limits" is set to')
    print(cas_labels)
    cassettes['cassette_label'] = cassettes.groupby('mO_name')['forward_start'].transform(get_cassette_label, cas_labels)

    # output invalid cassettes
    errors = cassettes['cassette_label'].isnull()
    if errors.any():
        print(cassettes[errors].to_string())
    else:
        print('All cassettes have a label')


    ############################################## SAVE ##################################################
    pickle_save(cassettes, cassettes_pickle_file)
    dataframe_out(cassettes, cassettes_text_file)
