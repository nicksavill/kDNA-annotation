import re
import numpy as np
import pathlib
from Bio.Seq import reverse_complement
from collections import Counter, OrderedDict

from .common import *

def identify_gRNAs(mini_align_file, maxi_align_file, minicircles, mRNAs, cassettes, filter):
    def extract_alignment(alignment, gRNAs, mRNAs, circle_name, minicircles, cas_dict, filter):
        mRNA_name, mRNA_align_o, gRNA_align_o, strand = alignment.split()
        if mRNA_name in mRNAs:
            mRNA_record = mRNAs[mRNA_name]
        else:
            return

        dd = {'c':'coding', 't':'template'}
        strand = dd[strand]

        # get the pairing for this alignment
        mRNA_seq = mRNA_align_o.replace('T', 'U')
        gRNA_seq = complement(gRNA_align_o)
        if gRNA_seq == None:
            return
        else:
            gRNA_seq = gRNA_seq.replace('T', 'U')
        pairing  = pairs(mRNA_seq, gRNA_seq)

        # remove maxicircle gene alignments 
        if ':' not in pairing:
            return

        # find 5'-most WC anchor region
        min_anchor_length = filter['min_anchor_length']
        min_gRNA_length = filter['min_gRNA_length']
        anchor_seq_regex = re.compile(filter['anchor_seq'])

        anchors = anchor_seq_regex.finditer(pairing)
        try:
            # get the longest 5'-most anchor
            anchor = [a.group(1) for a in anchors if len(a.group(1)) >= min_anchor_length][-1]
            anchor_pos = pairing.rfind(anchor)
            length = anchor_pos+len(anchor)
        except IndexError:
            return
        if length < min_gRNA_length:
            return

        # trim alignment so that it starts at the anchor
        mRNA_align = mRNA_align_o[:length]
        gRNA_align = gRNA_align_o[:length]
        pairing = pairing[:length]

        # find position on mRNA of gRNA and check edits occur there
        mRNA_start = mRNA_record['DNA_seq'].index(mRNA_align)
        if 't' not in mRNA_record['edits'][mRNA_start:mRNA_start+length]:
            return

        # find start of gRNA on minicircle and assign cassette_label or orphan or Maxi
        if circle_name.startswith('mO'):
            # find start of alignment on minicircle
            if strand == 'template':
                start = str(minicircles[circle_name].seq).index(gRNA_align)
            else:
                start = str(minicircles[circle_name].seq).index(reverse_complement(gRNA_align))
            # calculate end of gRNA on minicircle
            end = start+length
            # find cassettes this gRNA starts and ends in
            # cassette labels this gRNA starts and ends in
            c_start_label, c_end_label, cp = None, None, None
            if circle_name in cas_dict:
                for c in cas_dict[circle_name]:
                    if c[1] < start < c[4]:
                        c_start_label = c[0]
                        cp = c
                        break
                for c in cas_dict[circle_name]:
                    if c[1] < end < c[4]:
                        c_end_label = c[0]
                        break

                if c_start_label is not None and c_start_label is c_end_label:
                    # gRNA is fully inside a single cassette. Keep if start is in correct position
                    sp = filter['gRNA_search_region'][0]
                    ep = filter['gRNA_search_region'][1]
                    if strand == 'coding':
                        if cp[1]+sp <= start <= cp[1]+ep:
                            c_label = c_start_label
                            rel_start = start-cp[2] # position of gene relative to end of forward repeat
                        else:
                            return
                    else:
                        if cp[4]-ep <= end <= cp[4]-sp:
                            c_label = c_start_label
                            rel_start = cp[3]-end  # position of gene relative to start of reverse repeat
                        else:
                            return
                elif any((c_start_label, c_end_label)):
                    # gRNA starts or ends in cassette but crosses cassette boundary so ignore
                    return

                else:
                    # gRNA must be an orphan
                    if filter['allow_orphans']:
                        c_label = 'Orphan'
                        rel_start = 0
                    else:
                        return
            else:
                # no cassettes in mO which is probably truncated so ignore all alignments
                return

        else:
            c_label = 'Maxi'
            # set to zero so that filtering can work
            start = 0
            rel_start = 0

        # saved alignment information
        gRNA = OrderedDict()
        gRNA['mO_name']        = circle_name
        gRNA['cassette_label'] = c_label
        gRNA['strand']         = strand
        gRNA['length']         = length
        gRNA['rel_start']      = rel_start
        gRNA['circle_start']   = start
        gRNA['circle_end']     = start+length
        gRNA['mRNA_name']      = mRNA_name
        gRNA['product']        = mRNA_name.split('_')[0]
        gRNA['mRNA_start']     = mRNA_start
        gRNA['mRNA_end']       = mRNA_start+length
        gRNA['mRNA_seq']       = mRNA_record['edits'][mRNA_start:mRNA_start+length].replace('t', 'u').replace('T', 'U')
        gRNA['gRNA_seq']       = complement(gRNA_align).replace('T', 'U')
        gRNA['pairing']        = pairing
        gRNA['mismatches']     = gRNA['pairing'].count('.')
        gRNA['anchor_len']     = len(anchor_seq_regex.match(pairing[::-1]).group(1))
        gRNAs.append(gRNA)
 
    def drop_multi(gRNAs, min_multi_length):
        if len(gRNAs) == 1:
            return gRNAs
        gRNAs = gRNAs.sort_values(['length'], ascending=[False])
        last = 0 # last gRNA to use after sorted
        # if top two edit different versions of same mRNA keep both
        if ('_v' in gRNAs.iloc[0]['mRNA_name'] and '_v' in gRNAs.iloc[1]['mRNA_name'] and
            gRNAs.iloc[0]['product'] == gRNAs.iloc[0]['product']):
            last += 1
            # check here
            # last += 2
        # if next is a likely gRNA also keep
        if (last+1 < len(gRNAs) and gRNAs.iloc[last+1]['length'] >= min_multi_length):
            last += 1
        return gRNAs.iloc[0:last+1]

    def drop_comp(gRNAs):
        if len(gRNAs['strand'].unique()) == 1:
            return gRNAs
        # if gRNAs occur on both strands choose longest
        gRNAs = gRNAs.sort_values('length', ascending=False)
        last = 0 # last gRNA to use after sorted
        # if top two edit different versions of same mRNA keep both
        if ('_v' in gRNAs.iloc[0]['mRNA_name'] and '_v' in gRNAs.iloc[1]['mRNA_name'] and
            gRNAs.iloc[0]['product'] == gRNAs.iloc[0]['product']):
            last += 1
            # check here
            # last += 2
        return gRNAs.iloc[0:last+1]


    #######################################################################################################################
    # create a dictionary of cassettes for fast lookup
    cas_dict = {}
    for _, c in cassettes.iterrows():
        mO = c['mO_name']
        if mO not in cas_dict:
            cas_dict[mO] = []
        cas_dict[mO].append((c['cassette_label'], c['forward_start'], c['forward_end'], c['reverse_start'], c['reverse_end']))

    # read in the alignments from minicircle amd maxicircle alignments to edited mRNA
    # these files contain many duplicates, overlapping alignments, low quality false positives
    gRNAs = []
    with open(mini_align_file) as f:
        for name in f:
            extract_alignment(next(f), gRNAs, mRNAs, name[1:-1].split()[0], minicircles, cas_dict, filter)
    with open(maxi_align_file) as f:
        for name in f:
            extract_alignment(next(f), gRNAs, mRNAs, name[1:-1].split()[0], 'Maxicircle', cas_dict, filter)

    # convert alignments to dataframe for filtering
    gRNAs = pd.DataFrame(gRNAs)

    print(f'from alignment = {len(gRNAs)}')

    # filter out all duplicate alignments
    gRNAs.drop_duplicates(inplace=True)
    print(f'after drop duplicates = {len(gRNAs)}')

    # if >1 alignment begins or ends at same mapped position on an mRNA remove the smaller ones  
    # have to do this otherwise collapsing doesn't work if >1 gRNA in a cassette edits different parts of an mRNA
    gRNAs = gRNAs.groupby(['mO_name', 'cassette_label', 'strand', 'mRNA_name', 'circle_start', 'mRNA_end']).apply(drop_smaller).reset_index(drop=True)
    gRNAs = gRNAs.groupby(['mO_name', 'cassette_label', 'strand', 'mRNA_name', 'circle_start', 'mRNA_start']).apply(drop_smaller).reset_index(drop=True)
    gRNAs = gRNAs.groupby(['mO_name', 'cassette_label', 'strand', 'mRNA_name', 'circle_end', 'mRNA_start']).apply(drop_smaller).reset_index(drop=True)
    gRNAs = gRNAs.groupby(['mO_name', 'cassette_label', 'strand', 'mRNA_name', 'circle_end', 'mRNA_end']).apply(drop_smaller).reset_index(drop=True)
    print(f'after drop smaller sub-alignnments = {len(gRNAs)}')

    # trim the mismatches at the left-most end of the alignment
    gRNAs = gRNAs.apply(trim_mismatches, axis=1)

    # filter out low quality alignments
    mask1 = gRNAs['length'] < filter['min_gRNA_length']
    mask2 = gRNAs['anchor_len'] < filter['min_anchor_length']
    mask3 = gRNAs['mismatches'] > filter['max_mismatches']
    gRNAs = gRNAs[~(mask1 | mask2 | mask3)]
    # filter out low quality alignments, step 2, at least two of each
    # do this before multi-gRNAs otherwise we lose some gRNAs 
    print(f'after drop min quality gRNAs step 1 = {len(gRNAs)}')
    if filter['remove_low_quality']:
        mask1 = gRNAs['length'] <= filter['min_gRNA_length']+2
        mask2 = gRNAs['anchor_len'] <= filter['min_anchor_length']+1
        mask3 = gRNAs['mismatches'] >= max(0, filter['max_mismatches']-1)
        mask4 = gRNAs['mismatches'] >= filter['max_mismatches']
        mask5 = gRNAs['cassette_label'] == 'Orphan' 
        gRNAs = gRNAs[~((mask1 & mask2) | (mask1 & mask3) | (mask2 & mask3) | (mask4 & mask5))]
        print(f'after drop min quality gRNAs step 2 = {len(gRNAs)}')

    if filter['allow_orphans'] == 'auto':
        # remove orphan gRNAs not in expected position
        # get position of orphans
        orphan_positions = gRNAs.apply(orphan_position, args=(cassettes,), axis=1)
        mask1 = orphan_positions.isin(filter['orphan_positions'])
        mask2 = orphan_positions.isnull()
        gRNAs = gRNAs[mask1 | mask2]
        print(f'after drop incorrect orphans = {len(gRNAs)}')
    elif filter['allow_orphans'] == 'none':
        gRNAs = gRNAs[gRNAs['cassette_label'] != 'Orphan']
        print(f'after drop incorrect orphans = {len(gRNAs)}')

    # if two alignments exist in same cassette and one is much shorter then discard it
    gRNAs = gRNAs.groupby(['mO_name', 'cassette_label', 'strand']).apply(drop_multi, *(filter['min_multi_length'],)).reset_index(drop=True)
    print(f'after drop multi gRNAs = {len(gRNAs)}')

    # if alignments occur on both strands choose longest
    gRNAs = gRNAs.groupby(['mO_name', 'cassette_label']).apply(drop_comp).reset_index(drop=True)
    print(f'after drop complement gRNAs = {len(gRNAs)}')

    # remove other false positives that have not been removed at this point
    print(f'Removing identified false positives\n\t{filter["false_positives"]}')
    # boolean series of false positives to remove 
    y = gRNAs[['mO_name', 'cassette_label']].apply(lambda x: list(x) in filter['false_positives'], axis=1)
    gRNAs = gRNAs[~y]
    print(f'after drop false positives = {len(gRNAs)}')

    # collapse overlapping alignments into a longer alignment
    gRNAs = gRNAs.groupby(['mO_name', 'cassette_label', 'strand', 'mRNA_name']).apply(collapse, *(minicircles, mRNAs)).reset_index(drop=True)
    print(f'after collapsing = {len(gRNAs)}')

    # trim a second time as some collapsed alignments extend too far
    gRNAs = gRNAs.apply(trim_mismatches, axis=1)

    # assign the correct name to the alignment gRNA
    def assign_name(gRNA):
        strand_name = {'coding':'', 'template':'a'}
        if gRNA['mO_name'].startswith('Maxicircle'):
            return f"{gRNA['mO_name']}_g{gRNA['mRNA_name']}({gRNA['mRNA_start']+1}-{gRNA['mRNA_end']})"
        else:
            return f"{gRNA['mO_name']}({gRNA['cassette_label']}{strand_name[gRNA['strand']]})_g{gRNA['mRNA_name']}({gRNA['mRNA_start']+1}-{gRNA['mRNA_end']})"

    gRNAs['name'] = gRNAs.apply(assign_name, axis=1)
    return gRNAs

def identify_CSBs(minicircles, CSB_regexes):
    CSB1_regex  = re.compile(CSB_regexes['CSB1'])
    CSB2_regex  = re.compile(CSB_regexes['CSB2'])
    CSB3_regex  = re.compile(CSB_regexes['CSB3'])

    CSB1 = {}
    CSB2 = {}
    CSB3 = {}

    CSB1_count = Counter()
    CSB2_count = Counter()
    CSB3_count = Counter()
    for mO_name, minicircle_record in minicircles.items():
        seq = str(minicircle_record.seq)[:100]
        match = CSB1_regex.search(seq)
        if match:
            CSB1[mO_name] = {'start':match.start(0), 'end':match.end(0)}
            CSB1_count[match.group(0)] += 1

        match = CSB2_regex.search(seq)
        if match:
            CSB2[mO_name] = {'start':match.start(0), 'end':match.end(0)}
            CSB2_count[match.group(0)] += 1

        match = CSB3_regex.search(seq)
        if match:
            CSB3[mO_name] = {'start':match.start(0), 'end':match.end(0)}
            CSB3_count[match.group(0)] += 1

    print(f'CSB1: {len(CSB1)}, CSB2: {len(CSB2)}, CSB3: {len(CSB3)}')
    print(CSB1_count)
    print(CSB2_count)
    print(CSB3_count)
    return CSB1, CSB2, CSB3

def get_relative_position(gRNA, init_site):
    """ relative position of gRNA from 5' end of initiation site. For orphans distance is 0 """
    if gRNA['cassette_label'] in ['Orphan', 'Maxi']:
        return 0
    else:
        return gRNA['rel_start'] - init_site

def identify_gRNA_families(gRNAs, mRNAs, init_seq_len):
    """ assign gRNAs to gRNA families """
    gRNA_families = {'family_no':[], 'family_end':[], 'family_id':[]}
    strand_name = {'coding':'', 'template':'t'}
    index = []

    gRNAs['exp_mRNA_end'] = gRNAs['mRNA_end']+gRNAs['rel_pos']
    gRNAs['exp_mRNA_end'] = gRNAs['exp_mRNA_end'].astype('Int32')
    gRNAs['tmp'] = gRNAs.apply(lambda x: x['cassette_label']+strand_name[x['strand']], axis=1)

    for mRNA_name, mRNA in sorted(mRNAs.items()):
        # get all gRNAs with an for this mRNA
        g = gRNAs[gRNAs['mRNA_name'] == mRNA_name]

        # positions where the start of gRNAs align to mRNA + init_seq_len nt upstream of this position
        a = np.zeros(mRNA['length']+100)
        i = np.array(g['exp_mRNA_end']-1, dtype=int)
        for ii in range(init_seq_len):
            a[i-ii] = 1
        a = ''.join([str(int(i)) for i in a])
        g_end = 'exp_mRNA_end'

        tmp_g = []
        family_no = 0

        # find regions where groups of gRNAs anchor to same region of mRNA starting from 3' end of mRNA
        for m in re.finditer('1+', a):
            s, e = m.start(0), m.end(0)
            # the group of gRNAs encoded in the same cassette label and anchor to the same region of the mRNA
            anchor_group = g[(g[g_end] >= s) & (g[g_end] <= e)]

            if len(anchor_group) == 0:
                continue

            # for each cassette label of these gRNAs create a dictionary of cassette label and editing position
            cas_pos = {}
            for _, gRNA in anchor_group.iterrows():
                pos = gRNA['tmp']
                if pos not in cas_pos:
                    cas_pos[pos] = gRNA[g_end]
                cas_pos[pos] = max(gRNA[g_end], cas_pos[pos])

            # group gRNAs with the same cassette label ordered by editing position
            for pos, end in sorted(cas_pos.items(), key=lambda kv: kv[1]):
                group = anchor_group.query('tmp == @pos')
                index.extend(group.index.values)
                gRNA_families['family_no'].extend([family_no]*len(group))
                gRNA_families['family_end'].extend([end]*len(group))
                gRNA_families['family_id'].extend([f'{mRNA_name}-{pos}-{int(end)}']*len(group))
                tmp_g.append((family_no, end, f'{mRNA_name}-{pos}-{int(end)}'))
                family_no += 1

    gRNAs = gRNAs.drop(['tmp', 'exp_mRNA_end'], axis=1)
    return gRNAs.join(pd.DataFrame(gRNA_families, index=index))

def cassette_type(gRNAs, cassettes):
    """
        A cassette is canonical if canonical gRNA exists on either strand, otherwise non-canonical
    """
    # get all gRNAs and drop duplicates due to same gRNA editing different versions of same gene
    gRNAs = gRNAs[['mO_name', 'cassette_label']].drop_duplicates()
    gRNAs['exist'] = True

    l1 = len(cassettes)
    # merge gRNAs with cassettes to get cassette type, some cassettes may be duplicated
    cassettes = cassettes.merge(gRNAs, how='left')
    # assign cassette type based on existence of canonical gRNA 
    cassettes['type'] = cassettes['exist'].apply(lambda x: 'non-canonical' if x is np.nan else 'canonical')
    # cassettes = cassettes.drop(['strand'], axis=1)
    l2 = len(cassettes)
    assert(l1 == l2), f'{l1} {l2}'

    return cassettes


def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir, annotation_dir = get_directories(config)[1:3]

    minicircle_file       = f"{work_dir}/{config['minicircle clean fasta file']}"
    edited_mRNA_t_file    = f"{work_dir}/{config['edited mRNA with t fasta file']}"
    deletion_mRNA_file    = f"{work_dir}/{config['deletions mRNA text file']}"
    cassettes_pickle_file = f"{work_dir}/{config['cassettes pickle file']}"
    mini_align_file       = f"{work_dir}/{config['minicircle alignments file']}"
    maxi_align_file       = f"{work_dir}/{config['maxicircle alignments file']}"
    features_file         = f"{work_dir}/{config['features pickle file']}"
    hq_gRNAs_pickle_file  = f"{work_dir}/{config['high quality gRNAs pickle file']}"
    motifs_pickle_file    = f"{work_dir}/{config['motifs pickle file']}"
    gRNAs_text_file       = f"{annotation_dir}/{config['gRNAs text file']}"
    cassettes_text_file   = f"{annotation_dir}/{config['cassettes text file']}"

    pathlib.Path(annotation_dir).mkdir(parents=True, exist_ok=True) 


    ############################################### LOAD ###########################################
    minicircles  = get_minicircles(minicircle_file)
    mRNAs = get_mRNAs(edited_mRNA_t_file, deletion_mRNA_file)
    cassettes = pickle_load(cassettes_pickle_file)
    hq_gRNAs = pickle_load(hq_gRNAs_pickle_file)
    # load motif positions on each minicircle
    motif_positions = pickle_load(motifs_pickle_file)[0]


    ########################################## PARAMETERS #########################################
    # parameters for gRNA selection
    filter = config['all gRNAs filter']
    CSB_regexes = config['CSB regexes']
    # the position of the initiation site relative to the 18bp repeat 
    init_seq_len = config['initiation sequence length']

    # The region in which to search for gRNAs relative to the 5' end of the forward repeat
    # Use the HQ gRNA forward repeat positions  
    gRNA_search_start = config['upstream'] - motif_positions['forward repeat'].max()
    gRNA_search_end = gRNA_search_start + config['expected gRNA length']
    filter['gRNA_search_region'] = [gRNA_search_start, gRNA_search_end]
    print(f"gRNA search region relative to 5' end of forward repeat: {filter['gRNA_search_region']}")

    # The start of the initiation sequence uses the mode of the distances between the # 3' ends 
    # of the forward repeat motifs and the 5' ends of the initiation sequence motifs of HQ gRNAs
    # This is only used if expression information is not available 
    if config['have transcriptomics']:
        init_site = None
    else:
        init_site = (motif_positions['init sequence'] - motif_positions['forward repeat'] - config['repeat length']).mode().loc[0]
        print(f"Estimated 5' end of initiation site relative to 3' end of forward repeat: {init_site}")

    ##################### GET HQ ORPHAN POSITIONS FOR FILTERING ####################################
    if filter['allow_orphans'] == 'auto':
        # automcatically detect orphan position based on high quality gRNA orphans
        hq_gRNAs['orphan_position'] = hq_gRNAs.apply(orphan_position, args=(cassettes,), axis=1)
        orphans = hq_gRNAs[hq_gRNAs['orphan_position'].notnull()].drop(['mRNA_seq', 'gRNA_seq', 'sequence', 'product', 'mRNA_start', 'mRNA_end', 'anchor', 'mismatches'], axis=1)
        # add orphan positions to filter
        filter['orphan_positions'] = orphans['orphan_position'].unique()

        print('High quality orphans and their positions')
        print(orphans.to_string())
        print(filter['orphan_positions'])
    elif filter['allow_orphans'] == 'none':
        print('Orphans not allowed')


    ####################################### IDENTIFY CSBs AND gRNAS ###############################
    CSB1, CSB2, CSB3 = identify_CSBs(minicircles, CSB_regexes)

    gRNAs = identify_gRNAs(mini_align_file, maxi_align_file, minicircles, mRNAs, cassettes, filter)

    # Run the following if there are no transcriptomics data to identify expressed genes
    if not config['have transcriptomics']:
        # add anchors
        gRNAs, mRNAs = identify_anchors(gRNAs, mRNAs, filter)
        # relative postion of the gRNA to the initiation site 
        gRNAs['rel_pos'] = gRNAs.apply(get_relative_position, args=(init_site,), axis=1)
        # gRNA families 
        gRNAs = identify_gRNA_families(gRNAs, mRNAs, init_seq_len)
        # determine type of gRNA (canonical or non-canonical) in each cassette
        cassettes = cassette_type(gRNAs, cassettes)


    ############################################## SAVE ###########################################
    # save cassettes and gRNAs to annotation directory if there are no transcriptomics
    if not config['have transcriptomics']:
        dataframe_out(cassettes, cassettes_text_file)
        dataframe_out(gRNAs, gRNAs_text_file)

    # pickle everything together for annotation and analysis
    gzip_pickle_save([minicircles, mRNAs, CSB1, CSB2, CSB3, cassettes, gRNAs], features_file)
