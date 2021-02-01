import numpy as np
from scipy.stats import binom

from .common import *

def predict_expression(transcripts, init_site_range, p):
    """ Probability that at least the number of observed transcripts in the initiation site 
        is due to random chance assuming the cassette is not expressed
    """

    # def get_mode(df):
    #     """ get position of mode and number of transcripts in mode """
    #     common_pos = df['rel_pos'].mode()

    #     if len(common_pos) > 0:
    #         n = df[df['rel_pos'] == common_pos[0]]
    #         return pd.Series([common_pos[0], len(n)])
    #     else:
    #         return pd.Series([-1, 0])

    index = ['mO_name', 'cassette_label', 'strand']

    # keep only relevant information
    transcripts = transcripts[[*index, 'rel_pos']]
    grouped = transcripts.groupby(index)

    # set up masks for orphans/not orphans and transcripts in initiation site
    mask1 = transcripts['cassette_label'] != 'Orphan'
    mask2 = (init_site_range[0] <= transcripts['rel_pos']) & (transcripts['rel_pos'] <= init_site_range[1])

    # dataframe of cassette+strand with >0 transcripts 
    c = pd.DataFrame()
    # count total transcripts in each cassette on either strand
    c['transcripts_total'] = grouped['rel_pos'].count()

    # count transcripts in initiation site in each cassette on either strand
    # or return all transcripts for Orphans
    transcripts_init_site = transcripts[mask1 & mask2 | ~mask1]
    c['transcripts_init_site'] = transcripts_init_site.groupby(index)['rel_pos'].count()
    # set transcript count to 0 in cassettes without transcripts
    c['transcripts_init_site'].replace({np.nan:'0'}, inplace=True)
    c['transcripts_init_site'] = c['transcripts_init_site'].astype(int)

    # probability of observed number of transcripts starting in initiation site by random chance
    c['p-value'] = 1-binom.cdf(c['transcripts_init_site'], c['transcripts_total'], p)+binom.pmf(c['transcripts_init_site'], c['transcripts_total'], p)
    # predicted expression status bonferroni corrected by number of cassettes tested
    n = len(c.query('strand == "coding"'))
    c['expression'] = c['p-value'].apply(lambda x: 'expressed' if x < 0.05/n else 'non-expressed')

    # set expression to "outside" for any gRNAs whose initiation sequence is not in the initiation site
    # get position of mode and number of counts in mode (this is used for any gRNAs whose 
    # initiation sequence does not start in the initiation site)
    # c[['mode_pos', 'mode_count']] = grouped.apply(get_mode)
    # y = (c.index.get_level_values('strand') == 'coding') & (c['expression'] == "non-expressed") & ((c['mode_pos'] == 27) | (c['mode_pos'] == 33)) & (c['mode_count'] > 100)
    # c.loc[y, 'expression'] = 'outside' 

    return c[['transcripts_total', 'transcripts_init_site', 'p-value', 'expression']]

def expression(transcripts, init_site_range, init_seq_len, end_pos_percentile, p):
    def init_and_end(transcripts):
        # find the common position, the init_seq_len-nt sequence at that position and the number of transcripts at that position
        # find the initiation site position, the init_seq_len-nt sequence at that position and the number of transcripts at that position
        # predict the end of the sequenced gRNA by returning the end_pos_percentile 
        strand = transcripts.name[2]
        init = {'rel_start':0, 'common_start':0, 'init_seq':pd.NA, 'common_seq':pd.NA, 'transcripts_init_pos':0, 'transcripts_common_pos':0, 'rel_end':0}

        common_pos = transcripts['rel_pos'].mode()

        if len(common_pos) > 0:
            init['common_start'] = common_pos[0]
            mask = transcripts['rel_pos'] == common_pos[0]
            # the most common initiation sequence at the most common position
            if strand == 'coding':
                init['common_seq'] = transcripts[mask]['coding_seq'].str[:init_seq_len].mode()[0]
            else:
                init['common_seq'] = transcripts[mask]['coding_seq'].str[-init_seq_len:].mode()[0][::-1]
            init['transcripts_common_pos'] = mask.sum()

        # get transcription initiator position
        mask1 = transcripts['cassette_label'] != 'Orphan'
        mask2 = (init_site_range[0] <= transcripts['rel_pos']) & (transcripts['rel_pos'] <= init_site_range[1])
        common_pos = transcripts[mask1 & mask2 | ~mask1]['rel_pos'].mode()

        if len(common_pos) > 0:
            init['rel_start'] = common_pos[0]
            mask = transcripts['rel_pos'] == common_pos[0]
            # the most common initiation sequence at the most common position
            if strand == 'coding':
                init['init_seq'] = transcripts[mask]['coding_seq'].str[:init_seq_len].mode()[0]
            else:
                init['init_seq'] = transcripts[mask]['coding_seq'].str[-init_seq_len:].mode()[0][::-1]
            init['transcripts_init_pos'] = mask.sum()

        # a good predictor of the end of the expressed gRNA gene is the 90th percentile
        # end positions of transcripts with U tails.
        if strand == 'coding':
            transcripts = transcripts[transcripts['non_coding_3p_seq'].str[0] == 't']
        else:
            transcripts = transcripts[transcripts['non_coding_5p_seq'].str[-1] == 't']
        init['rel_end'] = transcripts['end_pos'].quantile(end_pos_percentile)
        return pd.Series(init)

    index = ['mO_name', 'cassette_label', 'strand']

    # predict expression in each cassette+strand
    expression_df = predict_expression(transcripts, init_site_range, p)

    # predict start and end of expressed gene and initiation sequence. include transcript counts
    grouped = transcripts.groupby(index)
    expression_df = expression_df.join(grouped.apply(init_and_end))
    expression_df['rel_end'] = expression_df['rel_end'].round().astype('Int64')

    # look for cassettes with high numbers of transcripts in the common_pos which aren't in the initiation site
    # expression_df = expression_df.reset_index()
    # mask1 = (expression_df['strand'] == 'coding') & (expression_df['expression'] == 'non-expressed')
    # mask2 = expression_df['common_start'] != expression_df['rel_start']
    # print(expression_df[mask1 & mask2].to_string())
    # exit()
    # modify the gRNAs whose initiation sequences are outside the initiation site
    # o = expression_df['expression'] == "outside"
    # expression_df.loc[o, 'rel_start'] = expression_df.loc[o, 'common_start']
    # expression_df.loc[o, 'init_seq'] = expression_df.loc[o, 'common_seq']
    # expression_df.loc[o, 'transcripts_init_pos'] = expression_df.loc[o, 'transcripts_common_pos']
    # expression_df.loc[o, 'expression'] = 'expressed'

    # remove unwanted columns
    expression_df = expression_df.drop(['common_start', 'common_seq', 'transcripts_init_pos', 'transcripts_common_pos'], axis=1)

    return expression_df

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir = get_directories(config)[1]

    transcripts_file       = f"{work_dir}/{config['transcripts pickle file']}"
    expression_pickle_file = f"{work_dir}/{config['expression pickle file']}"
    expression_text_file   = f"{work_dir}/{config['expression text file']}"


    ########################################## PARAMETERS #########################################
    p = config['p']
    init_seq_len = config['initiation sequence length']
    init_site_range = config['initiation site range']
    end_pos_percentile = config['end position percentile']/100


    ############################################### LOAD ###########################################
    transcripts = gzip_pickle_load(transcripts_file)


    ################################# DETERMINE EXPRESSION  ########################################
    # identify expressed genes and record number of transcripts
    genes = expression(transcripts, init_site_range, init_seq_len, end_pos_percentile, p)


    ##################################### SAVE ####################################################
    # output all cassettes+strands with transcripts
    # cassettes without any transcripts will not be in this file
    pickle_save(genes, expression_pickle_file)
    dataframe_out(genes, expression_text_file)
