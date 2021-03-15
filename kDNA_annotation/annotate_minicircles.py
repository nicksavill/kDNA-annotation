import pathlib
import pandas as pd
from collections import OrderedDict
from Bio.SeqFeature import SeqFeature, FeatureLocation

from .common import *

"""
Alignment files key:
Lines 1-4:  mRNA position (thousands, hundreds, tens, units) starting from 1
Line 5:     "A"s represent positions on the mRNA at which extender gRNAs anchour
            "I"s represent positions on the mRNA at which initiator gRNAs anchour
            "U"s represent positions on the mRNA at which extender gRNAs anchour but no 3' gRNA exists which edits the anchour 
Line 6:     "M"s represent positions of U-insertions not covered by gRNAs
Line 7:     "E"s represent positions on the mRNA covered by expressed gRNA editing regions
Line 8:     number of deletions to the right of the nucleotide below
Line 9:     edited mRNA sequence 5' to 3' (lowercase "u"s represent insertions)
Line 10:    protein sequence

For each gRNA:
Line 1:     name (mO_name(cassette position)_mRNA(start-end of alignment on mRNA)) 
            anchour represented by 
                "-" for extender gRNA
                ":" for initiator gRNA
                '.' for unanchoured gRNA
                '*' for undetermined extender or initiator
Line 2:     base-pairing
                "|":    Watson-Crick basepair
                ":":    GU basepair
                ".":    mismatch basepair
                "-":    gap (only in gapped dataset)
Line 3:     gRNA sequence 3' to 5'
Line 4:     extent of small RNA sequence

"""

def annotate(minicircles, CSB1, CSB2, CSB3, cassettes, gRNAs, genes, init_seq_len):
    cstrand = {'coding':1, 'template':-1}
    for mO_name, minicircle_record in minicircles.items():
        # minicircle_record.description = 'Trypanosoma brucei brucei strain AnTat1.1 90-13'
        minicircle_record.annotations['molecule_type'] = 'DNA'
        # minicircle_record.annotations['accession'] = ''
        # minicircle_record.annotations['version'] = ''
        # minicircle_record.annotations['keywords'] = ['kinetoplast DNA', 'kDNA', 'guide RNA', 'gRNA', '18bp inverted repeats']
        # minicircle_record.annotations['source'] = 'kinetoplast Trypanosoma brucei brucei'
        # minicircle_record.annotations['organism'] = 'Trypanosoma brucei brucei'
        # minicircle_record.annotations['topology'] = 'circular'
        # minicircle_record.annotations['date'] = str(datetime.datetime.now().strftime("%d-%b-%Y")).upper()
        # minicircle_record.annotations['taxonomy'] = ['Trypanosomatid', 'etc']
        # ref = Reference()
        # ref.authors = 'Cooper, S., Savill, N. J., Schnaufer, A.'
        # ref.title = 'kDNA genome T.b.brucei'
        # ref.journal = 'Unpublished'
        # minicircle_record.annotations['references'] = [ref]

        # as there can be multiple genes per cassette
        # we need to drop duplicate init_seq in the same cassette otherwise we can have
        # more than one initiation sequence per cassette
        def get_motif_start(sRNA):
            if sRNA['strand'] == 'coding':
                return sRNA['circle_start']
            else:
                return sRNA['circle_end']-init_seq_len

        if genes is not None:
            init_seq = genes.query('cassette_label != "Orphan"').dropna(subset=['circle_start']).copy()
            init_seq['start'] = init_seq.apply(get_motif_start, axis=1)
            init_seq = init_seq [['mO_name', 'start', 'strand']]
            init_seq = init_seq.drop_duplicates()
        else:
            init_seq = None

        features  = [SeqFeature(FeatureLocation(CSB1[mO_name]['start'], CSB1[mO_name]['end']), type='CSB1')]
        if mO_name in CSB2:
            features += [SeqFeature(FeatureLocation(CSB2[mO_name]['start'], CSB2[mO_name]['end']), type='CSB2')]
        features += [SeqFeature(FeatureLocation(CSB3[mO_name]['start'], CSB3[mO_name]['end']), type='CSB3')]

        features += [SeqFeature(FeatureLocation(i['forward_start'], i['reverse_end']),
            type='CASSETTE',
            qualifiers={'position':i['cassette_label']})
            for _, i in cassettes.query('mO_name == @mO_name').iterrows()]

        features += [SeqFeature(FeatureLocation(i['forward_start'], i['forward_end']),
            type='forward_repeat',
            # qualifiers={'seq':i['forward_seq']}
            )
            for _, i in cassettes.query('mO_name == @mO_name').iterrows()]

        if genes is not None:
            features += [SeqFeature(FeatureLocation(int(i['start']), int(i['start']+init_seq_len), strand=cstrand[i['strand']]),
                type='initiation',
                # qualifiers={'init_seq':i['init_seq']}
                )
                for _, i in init_seq.query('mO_name == @mO_name').iterrows()]

            features += [SeqFeature(FeatureLocation(int(i['circle_start']), int(i['circle_end']), strand=cstrand[i['strand']]),
                type='gene',
                qualifiers=OrderedDict([
                    ('mRNA_',             f"5'-{i['mRNA_seq']}-3'"),
                    ('align',             f"   {i['pairing']}   "),
                    ('seq__',             f"3'-{i['gRNA_seq']}-5'"),
                ])) for _, i in genes.query('mO_name == @mO_name and pairing == pairing').iterrows()]

        features += [SeqFeature(FeatureLocation(int(i['circle_start']), int(i['circle_end']), strand=cstrand[i['strand']]),
            type='gRNA',
            qualifiers=OrderedDict([
                ('name',              i['name']),
                # ('product', f"{i['mRNA_name']}"),
                # ('start', f"{i['mRNA_start']}"),
                # ('end', f"{i['mRNA_end']}"),
                # ('length', f"{i['length']}"),
                # ('anchor', f"{i['anchor']}"),
                # ('mismatches', f"{i['mismatches']}"),
                # ('pos', f"{i['cassette_label']}"),
                ('mRNA_',             f"5'-{i['mRNA_seq']}-3'"),
                ('align',             f"   {i['pairing']}   "),
                ('seq__',             f"3'-{i['gRNA_seq']}-5'"),
            ])) for _, i in gRNAs.query('mO_name == @mO_name').iterrows()]

        features += [SeqFeature(FeatureLocation(i['reverse_start'], i['reverse_end']),
            type='reverse_repeat',
            # qualifiers={'seq':i['reverse_seq']}
            )
            for _, i in cassettes.query('mO_name == @mO_name').iterrows()]

        minicircle_record.features = sorted(features, key=lambda x: x.location.start)

def output_edits(gRNAs, mRNAs, config, alignments_dir):

    a_type = {'extenderA':'_', 'extenderB':'_', 'initiator':'_', 'unanchored':'*'}

    if config['have transcriptomics']:
        gRNAs['gene_mRNA_end'] = gRNAs['mRNA_end']+gRNAs['rel_pos']
        gRNAs['gene_length'] = gRNAs['gene_rel_end']-gRNAs['gene_rel_start']
        gRNAs['gene_mRNA_start'] = gRNAs['gene_mRNA_end'] - gRNAs['gene_length']
        gRNAs['sort_pos'] = gRNAs['gene_mRNA_end']
    else:
        gRNAs['expression'] = 'unknown'
        gRNAs['sort_pos'] = gRNAs['mRNA_end']

    for mRNA_name, mRNA_record in mRNAs.items():
        mRNA_seq = mRNA_record['seq']

        # pack reads into rows
        row = 0
        alignments = [[]]
        # allow gRNAs to extend past the 3' end of a mRNA (to make small RNAs work)
        full_length = mRNA_record['length']+10

        # all gRNAs that edit this mRNA
        g = gRNAs.query('mRNA_name == @mRNA_name')

        nrows = 0
        rightmost = {}
        # go through groups from 3' to 5' end of mRNA 
        for family_no in sorted(g['family_no'].unique(), reverse=True):
            # all gRNAs in this family_no
            group_gRNAs = g[g['family_no'] == family_no].sort_values(['sort_pos', 'mRNA_end'], ascending=[False, False])

            # try putting gRNAs in topmost row first otherwise find row in which gRNA family_no will fit
            row = 0
            while True:
                for _, gRNA in group_gRNAs.iterrows():
                    # check if row exists, if not create it
                    if nrows == row:
                        nrows += 1
                        if nrows > 100:
                            print(f'Too many rows in edit alignment: {mRNA_name}')
                            exit()
                        rightmost[row] = full_length
                    row += 1
                    if gRNA['expression'] == 'expressed':
                        if gRNA['gene_mRNA_end'] >= rightmost[row-1]:
                            break
                    else:
                        if gRNA['mRNA_end'] >= rightmost[row-1]:
                            break
                else:
                    # all gRNAs fit so update alignments
                    break

            # add gRNAs in this group to rows
            row -= len(group_gRNAs)
            for _, gRNA in group_gRNAs.iterrows():
                if len(alignments) <= row:
                    alignments.append([])
                alignments[row].insert(0, gRNA)
                # get left most filled position in row so that the next group of gRNAs does not overlap
                # this includes the expressed small RNA 
                if gRNA['expression'] == 'expressed' and gRNA['gene_length'] is not pd.NA:
                    try:
                        rightmost[row] = min(gRNA['gene_mRNA_end']-gRNA['gene_length'], gRNA['mRNA_end']-gRNA['length'])
                    except TypeError:
                        print(gRNA)
                        exit()
                else:
                    rightmost[row] = gRNA['mRNA_end']-gRNA['length']
                row += 1

        if config['have transcriptomics']:
            x = ['-' for _ in range(full_length)]
            for row in alignments:
                for gRNA in row:
                    if gRNA['expression'] == 'expressed':
                        x[int(gRNA['gene_mRNA_end'])-1] = 'X'

            a = ['-' for _ in range(full_length)]
            for row in alignments:
                for gRNA in row:
                    if gRNA['expression'] == 'expressed':
                        a[int(gRNA['mRNA_end'])-1] = '@'

        else:
            a, x = [], []

        out = []
        for j in [1000, 100, 10, 1]:
            out.append(''.join([str((i//j)%10) if i%j == 0 else ' ' for i in range(1, len(mRNA_seq)+1)]))
        # out.append(''.join(['A' if i == 1 or i == 4 else 'I' if i == 2 else 'U' if i == 3 else '-' for i in mRNA_record['anchor']]))
        # out.append(''.join([str(int(i)) for i in mRNA_record['anchor_count']]))
        out.append(''.join(['M' if (k == 'u' or j != '-') and i == 0 else '-' for i, j, k in zip(mRNA_record['edited'], mRNA_record['deletions'], mRNA_seq)]))
        if config['have transcriptomics']:
            out.append(''.join(a))
            out.append(''.join(x))
        out.append(mRNA_record['deletions'])
        out.append(mRNA_seq)
        out.append(' '*mRNA_record['orf']+''.join([f'{i}  ' for i in mRNA_record['translate']]))
        for row in alignments:
            gRNA_name_align = [' ' for _ in range(full_length)]
            pairing_align   = [' ' for _ in range(full_length)]
            sequence_align  = [' ' for _ in range(full_length)]
            expression      = [' ' for _ in range(full_length)]
            for i in range(0, full_length, 10):
                gRNA_name_align[i] = '.'
            for gRNA in row:
                # find position of gRNA alignment to mRNA
                start = gRNA['mRNA_end']-gRNA['length']
                end   = gRNA['mRNA_end']

                # get anchor and its position
                # a = gRNA['anchor']
                a_end = end
                if gRNA['cassette_label'] == 'Maxi':
                    gRNA['cassette_label'] = ''

                # info = [gRNA['name'], str(gRNA['family_no']), f"{a_type[gRNA['anchor_type']]*gRNA['anchor_length']}"]
                info = []
                # if gRNA['expression'] == 'expressed':
                #     info += ['*']
                # info += [gRNA['name'], f"{a_type[gRNA['anchor_type']]*int(gRNA['anchor_length'])}"]
                # info += [gRNA['name']]
                info += [gRNA['name'], gRNA['family_id']]
                # info += [gRNA['name'], str(int(gRNA['family_no'])), str(gRNA['mRNA_end'])]
                # info += [gRNA['name'], str(gRNA['init_pos'])]
                gRNA_header = ' '.join(info)
                gRNA_name_align[a_end-len(gRNA_header):a_end] = list(gRNA_header)
                pairing_align[start:end]  = list(gRNA['pairing' ][:])
                sequence_align[start:end] = list(gRNA['gRNA_seq'][:])
                if gRNA['expression'] == 'expressed' and gRNA['gene_mRNA_start'] is not pd.NA:
                    try:
                        e_start = int(gRNA['gene_mRNA_start'])
                    except:
                        print(gRNA)
                        exit()
                    e_end = int(gRNA['gene_mRNA_end'])
                    expression[e_start:e_end] = ['-']*(e_end-e_start)

            out.append(''.join(gRNA_name_align))
            out.append(''.join(pairing_align))
            out.append(''.join(sequence_align))
            if config['have transcriptomics']:
                out.append(''.join(expression))
        with open(f'{alignments_dir}/{mRNA_name}.txt', 'w') as f:
            outs = '\n'.join(out)
            f.write(outs)

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
    config = load_config(config_file)
    work_dir, annotation_dir = get_directories(config)[1:3]

    if config['have transcriptomics']:
        features_file = f"{work_dir}/{config['features with expression pickle file']}"
    else:
        features_file = f"{work_dir}/{config['features pickle file']}"
    genbank_file = config['genbank text file']

    alignments_dir = f"{annotation_dir}/{config['alignments directory']}"
    genbank_dir = f"{annotation_dir}/{config['genbank directory']}"

    pathlib.Path(alignments_dir).mkdir(parents=True, exist_ok=True) 
    pathlib.Path(genbank_dir).mkdir(parents=True, exist_ok=True) 


    ########################################### PARAMETERS #########################################
    init_seq_len = config['initiation sequence length']

    ############################################### LOAD ###########################################
    if config['have transcriptomics']:
        minicircles, mRNAs, CSB1, CSB2, CSB3, cassettes, gRNAs, genes = gzip_pickle_load(features_file)
    else:
        minicircles, mRNAs, CSB1, CSB2, CSB3, cassettes, gRNAs = gzip_pickle_load(features_file)
        genes = None

    ##################### SAVE GENBANK FILE AND FULL ALIGNMENTS TO #################################
    if config['output genbank']:
        annotate(minicircles, CSB1, CSB2, CSB3, cassettes, gRNAs, genes, init_seq_len)
        minicircle_list = [v for _, v in sorted(minicircles.items())]
        SeqIO.write(minicircle_list, f'{genbank_dir}/{genbank_file}', 'genbank')

    output_edits(gRNAs, mRNAs, config, alignments_dir)

