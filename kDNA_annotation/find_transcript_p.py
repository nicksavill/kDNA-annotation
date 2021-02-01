import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

from .common import *

def estimate_normal_curve(transcripts, init_site_range, fit_range):
    """ probability of a randomly selected transcript being in the initiation site
        under the null hypothesis that the cassette is not expressed.
        Or, more precisely, that transcript positions are normally distributed
        in the cassette with mean mu and st.dev. sigma
    """

    # use cassette-associated transcripts outside of initiation site 
    # to minimise effect of expressed gRNAs on coding strand
    mask1 = (init_site_range[0] <= transcripts['rel_pos']) & (transcripts['rel_pos'] <= init_site_range[1])
    mask2 = transcripts['cassette_label'] != 'Orphan'
    transcripts = transcripts[~mask1 & mask2]

    # count number of transcripts at each position from forward repeat
    distribution = transcripts.groupby(['rel_pos'])['mO_name'].count()

    # discard ends as they are not normal, log counts and fit a quadratic
    # use positions 15 to 70 inclusive, minus init_length as not using initiation site
    init_length = init_site_range[1]- init_site_range[0] + 1
    ss = slice(fit_range[0], fit_range[1]-init_length)
    x = distribution.index[ss]
    y = np.log(distribution.values[ss])
    z = np.polyfit(x, y, 2)

    return distribution, x, y, z


def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir = get_directories(config)[1]

    transcripts_file              = f"{work_dir}/{config['transcripts pickle file']}"
    transcript_distribution_file  = f"{work_dir}/{config['transcript position distribution pickle file']}"


    ########################################## PARAMETERS #########################################
    init_site_range = config['initiation site range']
    fit_range = config['transcript position fit range']


    ############################################### LOAD ###########################################
    transcripts = gzip_pickle_load(transcripts_file)


    ######################################### DETERMINE p ##########################################
    distribution, x, y, z = estimate_normal_curve(transcripts, init_site_range, fit_range)

    # mean and standard dev of the fitted normal curve
    sigma = np.sqrt(-1/(2*z[0]))
    mu = -z[1]/(2*z[0])
    # probability a random transcript starts in the initiation site
    p = norm.cdf(init_site_range[1]+1, mu, sigma)-norm.cdf(init_site_range[0], mu, sigma)

    print(f'polynomial = {z}')
    print(f'mu = {mu}, sigma = {sigma}')
    print(f'probability that a transcript starts in the initiation site by random chance = {p}')


    # plot distribution of transcript relative positions to check for non-normal tails
    # change "transcript position fit range" in config file to adjust range of fit 
    f = np.poly1d(z)
    plt.bar(x, y, label='data')
    plt.plot(x, f(x), label='fit')
    plt.legend()
    plt.xlabel("Transcript position relative to 3' end of forward repeat")
    plt.ylabel('log frequency')
    plt.show()


    ##################################### SAVE #####################################################
    pickle_save(distribution, transcript_distribution_file)
