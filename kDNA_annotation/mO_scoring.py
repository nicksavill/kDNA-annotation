"""
Create scoring vectors for gRNAs and inverted repeats
Output minimum scores for identifying peaks in gRNA and repeat scores
from the hq gRNAs
"""

from operator import itemgetter
import numpy as np

from .common import *
from .savitzky_golay import savitzky_golay

def mO_nt_freq_scores(minicircles, nt_freqs, length):
    def get_score(mini_seq):
        missing_base = np.zeros(length)+np.log(0.001)
        return sum([nt_freqs.get(b, missing_base)[i] for (i, b) in enumerate(mini_seq)])/length

    scores = {'mO_name':[], 'score':[], 'smoothed':[], 'smoothed_gradient':[]}

    for mO_name, m in minicircles.items():
        print(mO_name)
        score = np.array([get_score(m.seq[i:i+length]) for i in range(len(m.seq)-length)])
        smoothed = savitzky_golay(score, 31, 5)
        gradient = savitzky_golay(score, 31, 5, 1)
        scores['mO_name'].append(mO_name)
        scores['score'].append(score)
        scores['smoothed'].append(smoothed)
        scores['smoothed_gradient'].append(gradient)
    return pd.DataFrame(scores)

def calculate_gRNA_nt_freqs(minicircles, gRNAs, feature, up, down):
    def get_seq(gRNA):
        s = gRNA[feature]
        return str(minicircles[gRNA['mO_name']].seq)[s-up:s+down]
    
    seqs = gRNAs.apply(get_seq, axis=1)

    nt_freq = pd.DataFrame(index=['A', 'C', 'G', 'T'])
    for i in range(up+down):
        c = seqs.str[i].value_counts(normalize=True)
        c.name = i
        nt_freq = nt_freq.join(c)

    # replace frequency of missing bases with 0.001 to prevent NaNs in frequency scores
    nt_freq = nt_freq.replace({np.nan:0.001})
    # log the frequencies
    log_nt_freq = np.log(nt_freq)
    return log_nt_freq.T.to_dict('list')

def get_gRNA_peaks(mO_gRNA_scores, num_cassettes, up):
    # determine the minimum score for probable peaks in gRNA scores
    gRNA_peaks = {'mO_name':[], 'score':[], 'position':[]}
    
    # get all peak scores for each 
    for _, mO in mO_gRNA_scores.iterrows():
        # identify the position and value of all peaks along a minicircle
        s = mO['smoothed']
        ds = mO['smoothed_gradient']
        x = np.arange(len(s))+up
        peaks = [(x1, s1) for (x1, s1, ds1, ds2) in zip(x[:-1], s[:-1], ds[:-1], ds[1:]) if ds1 > 0 and ds2 <= 0]

        for x, s in sorted(peaks, key=itemgetter(1), reverse=True):
            gRNA_peaks['mO_name'].append(mO['mO_name'])
            gRNA_peaks['score'].append(s)
            gRNA_peaks['position'].append(x)

    gRNA_peaks = pd.DataFrame(gRNA_peaks)
    scores = gRNA_peaks['score'].sort_values(ascending=False)
    min_gRNA_score = scores.iloc[num_cassettes]
    gRNA_peaks['quality'] = gRNA_peaks['score'] > min_gRNA_score
    return gRNA_peaks, min_gRNA_score

def get_score(seq, nt_freqs, length):
    missing_base = np.zeros(length)+np.log(0.001)
    return sum([nt_freqs.get(b, missing_base)[i] for (i, b) in enumerate(seq)])/length

def get_mO_repeat_scores(nt_freqs, minicircles, repeat_len):
    idx = []
    val = []
    for mO_name, mO in minicircles.items():
        # this minicircles sequence
        mO_seq = str(mO.seq)
        n = len(mO_seq)-repeat_len
        score = np.array([get_score(mO_seq[i:i+repeat_len], nt_freqs, repeat_len) for i in range(n)])
        idx.append(mO_name)
        val.append(score)

    return pd.Series(val, index=idx)

def convert_to_log(nt_freqs):
    # log nucleotide frequencies
    # if frequency is 0 set it to 0.001 to prevent -infinity
    for base, freqs in nt_freqs.items():
        nt_freqs[base] = [np.log(float(p)) if float(p) != 0.0 else np.log(0.001) for p in freqs]
    return nt_freqs

def get_min_score(mO_scores, num_cassettes):
    scores = []
    for score in mO_scores:
        scores.extend(list(score))
    return sorted(scores, reverse=True)[num_cassettes]

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir = get_directories(config)[1]

    minicircle_file      = f"{work_dir}/{config['minicircle clean fasta file']}"
    hq_gRNAs_pickle_file = f"{work_dir}/{config['high quality gRNAs pickle file']}"
    scores_file          = f"{work_dir}/{config['scores pickle file']}"
    motifs_pickle_file   = f"{work_dir}/{config['motifs pickle file']}"


    ########################################## PARAMETERS #########################################
    # number of nt upstream and downstream of start of gRNA to output for motif calling and nt frequency scoring 
    up = config['upstream']  
    down = config['downstream']

    # expected number of cassettes per minicircle. Used to get minimum score of gRNA peaks
    expected_cassettes = config['expected number of cassettes'] 

    repeat_len = config['repeat length']


    ####################################### LOAD SEQUENCES #########################################
    minicircles = get_minicircles(minicircle_file)
    hq_gRNAs = pickle_load(hq_gRNAs_pickle_file)

    # load nucleotide frequencies for forward and reverse repeat and their regular expressions
    motif_nt_freqs = pickle_load(motifs_pickle_file)[2]

    # the expected number of cassettes based on observation of histogram of 5' end s of HQ gRNAs
    num_cassettes = expected_cassettes * len(minicircles)


    ##################################### gRNA SCORES FOR EACH MINICIRCLE #############################
    # calculate gRNA nucleotide frequencies from 'up'nt upstream from start of hq gRNAs to 'down'nt downstream
    gRNA_frequencies = calculate_gRNA_nt_freqs(minicircles, hq_gRNAs, 'circle_start', up, down)

    # use these frequencies to produce scores along each minicircle to aid in cassette identification
    mO_gRNA_scores = mO_nt_freq_scores(minicircles, gRNA_frequencies, up+down)

    # get peaks from scores and the minimum threshold for the peaks
    gRNA_peaks, min_gRNA_score = get_gRNA_peaks(mO_gRNA_scores, num_cassettes, up)
    print(f'gRNA peak threshold {min_gRNA_score:.1f}')


    ##################################### REPEAT SCORES FOR EACH MINICIRCLE #############################
    # convert scores to log scale
    forward_nt_freq = convert_to_log(motif_nt_freqs['forward repeat'])
    reverse_nt_freq = convert_to_log(motif_nt_freqs['reverse repeat'])

    mO_forward_scores = get_mO_repeat_scores(forward_nt_freq, minicircles, repeat_len)
    mO_reverse_scores = get_mO_repeat_scores(reverse_nt_freq, minicircles, repeat_len)

    # determine the minimum scores for probable peaks in repeat scores
    min_forward_score = get_min_score(mO_forward_scores, num_cassettes)
    min_reverse_score = get_min_score(mO_reverse_scores, num_cassettes)
    print(f'forward repeat peak threshold {min_forward_score:.1f}')
    print(f'reverse repeat peak threshold {min_reverse_score:.1f}')


    ####################################### SAVE ####################################################
    # save for later use in identify_cassettes.py
    i = [
        mO_gRNA_scores,
        mO_forward_scores,
        mO_reverse_scores,
        gRNA_peaks,
        min_gRNA_score,
        min_forward_score,
        min_reverse_score
    ]

    pickle_save(i, scores_file)

