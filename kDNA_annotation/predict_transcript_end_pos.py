import matplotlib.pyplot as plt

from .common import *

def predict_end_position(transcripts, rel_end_statistics):
    """
    """

    def init_and_end(transcripts, rel_end_statistics):
        strand = transcripts.name[2]
        init = {'mean_tail':0, 'mode_tail':0, '50%_tail':0, '90%_tail':0, '95%_tail':0, '75%_tail':0,
                'mean_all':0, 'mode_all':0, '50%_all':0, '90%_all':0, '95%_all':0, '75%_all':0}

        init = dict([(s, 0) for s in rel_end_statistics])

        init['75%_all'] = transcripts['end_pos'].quantile(0.75)
        init['90%_all'] = transcripts['end_pos'].quantile(0.90)
        init['95%_all'] = transcripts['end_pos'].quantile(0.95)
        init['mean_all'] = transcripts['end_pos'].mean()
        init['50%_all'] = transcripts['end_pos'].median()
        cp = transcripts['end_pos'].mode()
        if len(cp) > 0:
            init['mode_all'] = cp[0]
        else:
            init['mode_all'] =  np.nan

        # use transcripts with poly-U tail
        if strand == 'coding':
            transcripts = transcripts[transcripts['non_coding_3p_seq'].str[0] == 't']
        else:
            transcripts = transcripts[transcripts['non_coding_5p_seq'].str[-1] == 't']
        init['75%_tail'] = transcripts['end_pos'].quantile(0.75)
        init['90%_tail'] = transcripts['end_pos'].quantile(0.90)
        init['95%_tail'] = transcripts['end_pos'].quantile(0.95)
        init['50%_tail'] = transcripts['end_pos'].median()
        init['mean_tail'] = transcripts['end_pos'].mean()
        cp = transcripts['end_pos'].mode()
        if len(cp) > 0:
            init['mode_tail'] = cp[0]
        else:
            init['mode_tail'] =  np.nan

        return pd.Series(init)

    index = ['mO_name', 'cassette_label', 'strand']

    transcripts = transcripts[index+['end_pos', 'non_coding_3p_seq', 'non_coding_5p_seq']]
    return transcripts.groupby(index).apply(init_and_end, rel_end_statistics)


def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir = get_directories(config)[1]

    transcripts_file = f"{work_dir}/{config['transcripts pickle file']}"
    features_file    = f"{work_dir}/{config['features pickle file']}"


    ############################################### LOAD ###########################################
    transcripts = gzip_pickle_load(transcripts_file)
    gRNAs = gzip_pickle_load(features_file)[-1]


    #################### COMPARE END POSITION STATS ################################################
    rel_end_statistics = ['mean_tail', 'mode_tail', '50%_tail', '75%_tail', '90%_tail', '95%_tail', 
        'mean_all', 'mode_all', '50%_all', '75%_all', '90%_all', '95%_all']

    end_pos = predict_end_position(transcripts, rel_end_statistics)

    # get the column names of the different methods of obtaining end position
    labels = ['Mean', 'Mode', '50th percentile', '75th percentile', '90th percentile', '95th percentile']

    # select only expressed gRNAs and merge with end rel_end_statistics
    end_pos = end_pos.reset_index().merge(gRNAs, how='outer')

    stdev = {}

    fig, axes = plt.subplots(2, 6, sharex=True, sharey=True, figsize=(20, 10))
    for p, axis in zip(rel_end_statistics, axes.flatten()):
        # distance between position statistic of 3' end of gene and 3' end of aligned gRNA 
        x = end_pos[p] - (end_pos['rel_start'] + end_pos['length'])
        axis.hist(x, bins=range(-30, 30), align='left')
        axis.axvline(0, ls='--')
        axis.set_title(p)
        axis.annotate(f'{x.describe().iloc[1:].round(1).to_string()}', (0, 150), fontsize=10, fontfamily='monospace', va='top', ha='center')
        stdev[p] = x.std()
    axes[0, 0].set_ylabel('Number of expressed\ncanonical gRNAs')
    axes[1, 0].set_ylabel('Number of expressed\ncanonical gRNAs')
    for i in range(6):
        axes[1, i].set_xlabel(f'{labels[i]}\ndistance (nt)')
        axes[0, i].set_title('With U-tail')
        axes[1, i].set_title('Without U-tail')
    fig.suptitle("Comparison of statistics for defining the 3' end of expressed gRNA genes")

    stdev = pd.Series(stdev)
    idx = stdev.idxmin()
    print(f'Minimum variance statistic is {idx}: {stdev[idx]:.2f}')
    print(f'set "end position percentile" to {idx[:2]}')

    plt.show()
    return fig