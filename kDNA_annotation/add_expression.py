"""
    Add expression information to gRNAs and cassettes
"""

from copy import copy
from operator import itemgetter
import re
import numpy as np

from .common import *

def canonical_gRNA_expression(gRNAs, cassettes, expression, filter):
    def get_rel_pos(gRNA, min_anchor_length):
        """ relative position of gRNA from 5' end of initiation sequence
            For orphans distance is already given by minus rel_start
        """
        if gRNA['cassette_label'] == 'Orphan':
            return -gRNA['rel_start']

        if gRNA['strand'] == 'coding':
            rel_pos = gRNA['circle_start']-gRNA['forward_end']-gRNA['gene_rel_start']
        else:
            rel_pos = gRNA['reverse_start']-gRNA['gene_rel_start']-gRNA['circle_end']
        if rel_pos is pd.NA:
            rel_pos = 0

        if rel_pos < 0:
            # find position of first non-WC in pairing
            # If the resulting shortening of the alignment causes the anchor to be
            # less than the min anchor length, make rel_pos just past the non-WC bp
            match = mm_regex.search(gRNA['pairing'][::-1])
            mm_dist = match.start(0)
            if mm_dist + rel_pos < min_anchor_length:
                rel_pos = -(mm_dist+1)
        return rel_pos

    def trim_gRNA(gRNA):
        # adjust gRNA values for new rel_pos
        d = int(gRNA['rel_pos'])
        trimmed = {}

        if gRNA['strand'] == 'coding':
            trimmed['circle_start'] = gRNA['circle_start'] - d
            trimmed['circle_end'] = gRNA['circle_end']
        else:
            trimmed['circle_start'] = gRNA['circle_start']
            trimmed['circle_end'] = gRNA['circle_end'] + d
        trimmed['pairing'] = gRNA['pairing'][:d]
        trimmed['gRNA_seq'] = gRNA['gRNA_seq'][:d]
        trimmed['mRNA_seq'] = gRNA['mRNA_seq'][:d]
        trimmed['length'] = gRNA['length'] + d
        trimmed['mRNA_end'] = gRNA['mRNA_end'] + d
        return pd.Series(trimmed)

    index = ['pairing', 'gRNA_seq', 'mRNA_seq', 'length', 'mRNA_end', 'circle_start', 'circle_end']
    mm_regex = re.compile(r'[\.:]')

    # merge gRNAs with expression data to get transcription information
    gRNAs = gRNAs.merge(expression, on=['mO_name', 'cassette_label', 'strand'], how='left', suffixes=(None, '_g'))

    # rename columns of the expressed gene relative start and end
    gRNAs = gRNAs.rename(columns={'rel_start_g':'gene_rel_start', 'rel_end':'gene_rel_end'})
    # convert some columns to nullable integer

    # set transcript values to 0 if there are no transcripts 
    gRNAs['transcripts_total'] = gRNAs['transcripts_total'].replace({np.nan:0})
    gRNAs['transcripts_total'] = gRNAs['transcripts_total'].astype(int)
    gRNAs['transcripts_init_site'] = gRNAs['transcripts_init_site'].replace({np.nan:0})
    gRNAs['transcripts_init_site'] = gRNAs['transcripts_init_site'].astype(int)
    gRNAs['expression'] = gRNAs['expression'].replace({np.nan:'non-expressed'})
    gRNAs['p-value'] = gRNAs['p-value'].replace({np.nan:1})
    # calculate rel_pos from start of anchor to initiation position
    # if it is negative, adjust gRNA so that it starts at or past the initiation position
    # keeping the anchor at the start of the gene
    # don't trim Orphans
    # merge gRNAs with cassettes to get forward_end and reverse_start to calculate rel_pos of anchor to init_pos
    gRNAs = gRNAs.merge(cassettes[['mO_name', 'cassette_label', 'forward_end', 'reverse_start']], how='left')

    for c in ['forward_end', 'reverse_start', 'gene_rel_start', 'transcripts_total', 'transcripts_init_site']:
        gRNAs[c] = gRNAs[c].round().astype('Int64')

    gRNAs['rel_pos'] = gRNAs.apply(get_rel_pos, axis=1, args=(filter['min_anchor_length'],))

    if filter['trim_to_init']:
        # trim gRNAs with negative rel_pos
        # drop columns that will be changed and then join the updated columns
        c = (gRNAs['rel_pos'] < 0) & (gRNAs['cassette_label'] != 'Orphan')
        x = gRNAs.loc[c]
        # only do if len(x) > 0 otherwise an error occurs
        if len(x):
            gRNAs.loc[c] = x.drop(index, axis=1).join(x.apply(trim_gRNA, axis=1))
        # updated rel_pos
        gRNAs['rel_pos'] = gRNAs.apply(get_rel_pos, axis=1, args=(filter['min_anchor_length'],))
    
    # drop unnecessary columns
    return gRNAs.drop(['forward_end', 'reverse_start'], axis=1)

def get_expressed_genes(gRNAs, cassettes, mRNAs, minicircles, expression):
    """
    expressed gRNAs start at the initiation position and end at the 90th-percentile
    position of transcripts with a u-tail. 
    If no transcripts have a u-tail there will no sequences and pairing.
    On coding strand canonical and non-canonical can be expressed gRNAs.
    On template strand only canonical can be expressed gRNAs.
    """
    def alignment(gene, gRNA, d, da, gap=0):
        def score(gene, d, offset=0):
            # score is the number of WC from position d until 5 nt downstream from init seq
            WC_regex = re.compile(r'\|*')
            p = gene['pairing'][gene['length']-d+1:gene['length']-5+offset]
            match = WC_regex.match(p)
            if match:
                return len(match.group(0))
                # return p, match.group(0), len(match.group(0))
            else:
                return 0

        # copy expressed gRNA to make a new copy for saving
        gene = copy(gene)
        gene['mRNA_end_x'] = int(gRNA['mRNA_end'] + da + gap)
        gene['mRNA_start'] = gene['mRNA_end_x'] - gene['length'] - gap
        gene['mRNA_seq']   = mRNAs[gRNA['mRNA_name']]['seq'][gene['mRNA_start']:gene['mRNA_end_x']]

        if gRNA['strand'] == 'coding':
            gene['gRNA_seq'] = str(minicircles[gRNA['mO_name']].seq)[gene['circle_start_x']:gene['circle_end']].replace('T', 'U')[::-1]
        else:
            gene['gRNA_seq'] = complement(str(minicircles[gRNA['mO_name']].seq)[gene['circle_start_x']:gene['circle_end']]).replace('T', 'U')

        if gap == 1:
            # add gap into gRNA sequence
            gene['gRNA_seq'] = gene['gRNA_seq'][:-d] + '-' + gene['gRNA_seq'][-d:]
        elif gap == -1:
            # add gap into mRNA sequence
            gene['mRNA_seq'] = gene['mRNA_seq'][:-d+1] + '-' + gene['mRNA_seq'][-d+1:]

        # get pairing alignment
        gene['pairing'] = pairs(gene['mRNA_seq'].upper(), gene['gRNA_seq'])
        return gene, score(gene, d)

    def get_gene_alignment(gene):
        gene_mod = {}
        # if gene['circle_start'] is pd.NA:
        #     return pd.Series(gene_mod)

        if gene['cassette_label'] == 'Orphan':
            try:
                gene_mod['circle_start_x'] = gene['circle_start']+gene['rel_start']
            except:
                print(gene)
                exit()
        else:
            if gene['strand'] == 'coding':
                gene_mod['circle_start_x'] = gene['forward_end']+gene['rel_start']
            else:
                gene_mod['circle_end'] = gene['reverse_start']-gene['rel_start']

        if gene['rel_end'] is not pd.NA:
            if gene['cassette_label'] == 'Orphan':
                gene_mod['circle_end'] = gene['circle_start']+gene['rel_end']
            else:
                if gene['strand'] == 'coding':
                    gene_mod['circle_end'] = gene['forward_end']+gene['rel_end']
                else:
                    gene_mod['circle_start_x'] = gene['reverse_start']-gene['rel_end']
            gene_mod['length'] = abs(gene_mod['circle_end']-gene_mod['circle_start_x'])

            # if canonical gene exists than get new pairing
            if isinstance(gene['mRNA_name'], str):
                # distance to anchor from initiation sequence
                da = gene['rel_pos']

                # nogap default alignment
                nogap_gene, nogap_score = alignment(gene_mod, gene, da, da) 
                possible_gRNAs = []

                if da > 7:
                    dmm = len(nogap_gene['pairing'][-da:])-nogap_gene['pairing'][-da:].find('.')
                    for d in {da, dmm}:
                        # deletion in gene: add gap just before anchor in gene
                        possible_gRNAs.append(alignment(gene_mod, gene, d, da, gap=1))

                        # insertion in gene: add gap just before anchor in mRNA
                        possible_gRNAs.append(alignment(gene_mod, gene, d, da, gap=-1))

                # only keep genes with a score > nogap_score+1
                better_gRNAs = [i for i in possible_gRNAs if i[1] > nogap_score+2]
                if len(better_gRNAs) == 0:
                    gene_mod = nogap_gene
                else:
                    gene_mod = max(better_gRNAs, key=itemgetter(1))[0]
            else:
                # get gene_mod sequence of non-canonical expressed gene_mod
                gene_mod['gRNA_seq'] = str(minicircles[gene['mO_name']].seq)[gene_mod['circle_start_x']:gene_mod['circle_end']].replace('T', 'U')[::-1]
        return pd.Series(gene_mod)

    # filter out non-expressed gRNAs
    genes = expression.query('expression == "expressed"')
    # merge with canonical gRNAs to get gRNA info
    genes = genes.merge(gRNAs[['mO_name', 'cassette_label', 'strand', 'rel_pos', 'mRNA_name', 'mRNA_end', 'circle_start']], on=['mO_name', 'cassette_label', 'strand'], how='left')
    # remove expressed non-canonical on template strand 
    genes = genes[~((genes['strand'] == 'template') & (genes['mRNA_name'].isnull()))]
    # merge genes with cassettes to get forward_end and reverse_start to calculate distance of anchor to rel_start
    genes = genes.merge(cassettes[['mO_name', 'cassette_label', 'forward_end', 'reverse_start']], how='left')
    # convert to nullable integer
    for c in ['rel_pos', 'mRNA_end', 'circle_start', 'forward_end', 'reverse_start']:
        genes[c] = genes[c].round().astype('Int64')

    # get new expressed gRNA sequences, pairings and length
    genes = genes.join(genes.apply(get_gene_alignment, axis=1))
    # drop unnecessary columns
    genes = genes.drop(['forward_end', 'reverse_start', 'rel_pos', 'expression', 'mRNA_end', 'circle_start'], axis=1)
    # rename mRNA_end_x to mRNA_end
    genes = genes.rename(columns={'mRNA_end_x':'mRNA_end', 'circle_start_x':'circle_start'})

    # set length and type of gene
    genes['length'] = genes['gRNA_seq'].str.len()
    genes['type'] = genes['mRNA_name'].apply(lambda x: 'canonical' if isinstance(x, str) else 'non-canonical')
    for c in [ 'circle_start','circle_end', 'length', 'mRNA_start', 'mRNA_end']:
        genes[c] = genes[c].round().astype('Int64')
    return genes

def cassette_type_expression(gRNAs, cassettes, expression):
    """
        A cassette is canonical if canonical gRNA exists on either strand, otherwise non-canonical
        A cassette is expressed if coding strand expressed or 
        canonical gRNA on either strand is expressed, otherwise nonexpressed
    """
    # get expression of all gRNAs and drop duplicates due to same gRNA editing different versions of same gene
    gRNAs = gRNAs[['mO_name', 'cassette_label', 'expression']].drop_duplicates()

    l1 = len(cassettes)
    # merge gRNAs with cassettes to get cassette type, some cassettes may be duplicated
    # only cassettes with canonical gRNAs will have expression data at this point
    cassettes = cassettes.merge(gRNAs, how='left')
    # assign cassette type based on existence of canonical gRNA 
    cassettes['type'] = cassettes['expression'].apply(lambda x: 'non-canonical' if x is np.nan else 'canonical')
    # cassettes = cassettes.drop(['strand'], axis=1)
    l2 = len(cassettes)
    assert(l1 == l2), f'{l1} {l2}'

    # select non-canonical cassettes
    mask = cassettes['type'] == 'non-canonical'
    # get mO name and cassette pos of noncanonical cassettes
    nonc = cassettes[mask][['mO_name', 'cassette_label']]
    idx = nonc.index

    # only use expression information for coding strand in cassette as there are too many false positive "expressed" template strand cassettes
    expression = expression.reset_index().query('strand == "coding"')[['mO_name', 'cassette_label', 'expression']]
    nonc = nonc.merge(expression, how='left').set_index(idx)
    cassettes.loc[mask, ['expression']] = nonc[['expression']]

    # 8 cassettes have no transcripts (so not in expression dataframe) so call them nonexpressed
    cassettes['expression'].replace({np.nan:'non-expressed'}, inplace=True)
    l3 = len(cassettes)
    assert(l2 == l3)
    return cassettes

def identify_expressed_gRNA_families(gRNAs, mRNAs, init_seq_len):
    """ assign gRNAs to editing groups based on anchor domains and cassette position """
    gRNA_families = {'family_no':[], 'family_end':[], 'family_id':[]}
    strand_name = {'coding':'', 'template':'t'}
    index = []

    gRNAs['gene_mRNA_end'] = gRNAs['mRNA_end']+gRNAs['rel_pos']
    gRNAs['gene_mRNA_end'] = gRNAs['gene_mRNA_end'].astype('Int32')
    gRNAs['tmp'] = gRNAs.apply(lambda x: x['cassette_label']+strand_name[x['strand']], axis=1)

    for mRNA_name, mRNA in sorted(mRNAs.items()):
        # get all gRNAs with an init_pos for this mRNA
        # nonexpressed gRNAs can be in an editing group if they have a init_seq. this is because
        # they have transcripts in the init_position but not enough to be called expressed
        # gRNAs without an init_seq have no transcripts within the initiation site
        # these are added to a group below
        mask1 = gRNAs['mRNA_name'] == mRNA_name
        mask2 = gRNAs['init_seq'].notnull()
        g = gRNAs[mask1 & mask2]

        # positions where the start of expressed gRNAs align to mRNA
        a = np.zeros(mRNA['length']+100)
        i = np.array(g['gene_mRNA_end']-1, dtype=int)
        for ii in range(init_seq_len):
            a[i-ii] = 1
        a = ''.join([str(int(i)) for i in a])
        g_end = 'gene_mRNA_end'

        tmp_g = []
        family_no = 0

        # find regions where groups of gRNAs anchor to mRNA starting from 3' end of edited mRNA
        for m in re.finditer('1+', a):
            s, e = m.start(0), m.end(0)
            # s, e = m.start(0)+1, m.end(0)-1
            # get all gRNAs that anchor at this region
            anchor_group = g[(g[g_end] >= s) & (g[g_end] <= e)]

            if len(anchor_group) == 0:
                continue

            # for each cassette position of these gRNAs create a dictionary of cassette position and editing position
            cas_pos = {}
            for _, gRNA in anchor_group.iterrows():
                pos = gRNA['tmp']
                if pos not in cas_pos:
                    cas_pos[pos] = gRNA[g_end]
                cas_pos[pos] = max(gRNA[g_end], cas_pos[pos])

            # group gRNAs with the same cassette position ordered by editing position
            for pos, end in sorted(cas_pos.items(), key=lambda kv: kv[1]):
                group = anchor_group.query('tmp == @pos')
                index.extend(group.index.values)
                gRNA_families['family_no'].extend([family_no]*len(group))
                gRNA_families['family_end'].extend([end]*len(group))
                gRNA_families['family_id'].extend([f'{mRNA_name}-{pos}-{int(end)}']*len(group))
                tmp_g.append((family_no, end, f'{mRNA_name}-{pos}-{int(end)}'))
                family_no += 1

        # gRNAs without an init_seq
        mask2 = gRNAs['init_seq'].isnull()
        unknown = gRNAs[mask1 & mask2]
        # for each unknown gRNA
        for idx, gRNA in unknown.iterrows():
            # search for a group that ends just after mRNA_end of this unknown gRMA
            for f_no, gene_mRNA_end, family_id in sorted(tmp_g, key=itemgetter(1)):
                [g_mRNA_name, g_pos, g_end] = family_id.split('-')
                if g_mRNA_name == mRNA_name and gRNA['mRNA_end']-1 <= gene_mRNA_end and gRNA['cassette_label'] == g_pos:
                    index.append(idx)
                    gRNA_families['family_no'].append(f_no)
                    gRNA_families['family_end'].append(gene_mRNA_end)
                    gRNA_families['family_id'].append(f'{family_id}')
                    break
            else:
                # no suitable gRNA found, so make a unique family for this non-expressed gRNA 
                index.append(idx)
                gRNA_families['family_no'].append(family_no)
                gRNA_families['family_end'].append(gRNA['mRNA_end'])
                gRNA_families['family_id'].append(f'{mRNA_name}-{gRNA["cassette_label"]}-{gRNA["mRNA_end"]}')
                family_no += 1

    gRNAs = gRNAs.drop(['tmp'], axis=1)
    gRNAs = gRNAs.join(pd.DataFrame(gRNA_families, index=index))
    gRNAs['family_no'] = gRNAs['family_no'].astype('Int64')
    gRNAs['family_end'] = gRNAs['family_end'].astype('Int64')
    return gRNAs

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir, annotation_dir = get_directories(config)[1:3]

    expression_pickle_file = f"{work_dir}/{config['expression pickle file']}"
    features_file          = f"{work_dir}/{config['features pickle file']}"
    features_with_exp_file = f"{work_dir}/{config['features with expression pickle file']}"
    genes_text_file        = f"{annotation_dir}/{config['expressed genes text file']}"
    gRNAs_text_file        = f"{annotation_dir}/{config['gRNAs with expression text file']}"
    cassettes_text_file    = f"{annotation_dir}/{config['cassettes with expression text file']}"


    ########################################## PARAMETERS #########################################
    filter = config['all gRNAs filter']
    init_seq_len = config['initiation sequence length']
    if not config['have transcriptomics']:
        print('The parameters "have transcriptomics" is set to "no" in config file.')
        print('add_expression.py can only be run if this parameter is "yes".')
        print('Re-run pipeline from identify_CSB_gRNAs.py" onward if transcriptomics are available.')
        exit()

    ############################################### LOAD ###########################################
    minicircles, mRNAs, CSB1, CSB2, CSB3, cassettes, gRNAs = gzip_pickle_load(features_file)
    expression = pickle_load(expression_pickle_file)

    if 'anchor_type' in gRNAs.columns:
        print(f'The file {features_file} has been created with "have transcriptomics" set to "no"')
        print('Set it to "yes" and re-run pipeline from identify_CSB_gRNAs.py" onward if transcriptomics is available.')
        exit()

    ###################### ADD EXPRESSION DATA TO gRNAs AND CASSETTES ###########################
    gRNAs = canonical_gRNA_expression(gRNAs, cassettes, expression, filter)
    # identify anchor after expression in case 5' ends need to be trimmed to start of initiation sequence
    gRNAs, mRNAs = identify_anchors(gRNAs, mRNAs, filter)
    gRNAs = identify_expressed_gRNA_families(gRNAs, mRNAs, init_seq_len)
    genes = get_expressed_genes(gRNAs, cassettes, mRNAs, minicircles, expression)
    cassettes = cassette_type_expression(gRNAs, cassettes, expression)


    ##################################### SAVE ####################################################
    dataframe_out(genes, genes_text_file)
    dataframe_out(gRNAs, gRNAs_text_file)
    dataframe_out(cassettes, cassettes_text_file)

    # save for annotation and analysis
    gzip_pickle_save([minicircles, mRNAs, CSB1, CSB2, CSB3, cassettes, gRNAs, genes], features_with_exp_file)
