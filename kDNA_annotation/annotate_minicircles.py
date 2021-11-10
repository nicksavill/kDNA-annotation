import pathlib
import pandas as pd
import datetime
from collections import OrderedDict
from Bio.SeqFeature import SeqFeature, FeatureLocation, Reference
from operator import itemgetter
from os import system

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

def submission(minicircles, CSB1, CSB2, CSB3, cassettes, gRNAs, genes, init_seq_len, genbank_dir):
    def make_locus_tag(gene):
        return f'GT92_{gene["mO_name"]}_{gene["cassette_label"]}_{gene["mRNA_name"]}'

    features = dict([(i, []) for i in minicircles])
    
    # extract all canonical gRNA info
    gRNAs = gRNAs[gRNAs['mO_name'] != 'Maxicircle']
    genes['expression'] = 'expressed'
    genes['gene_circle_start'] = genes['circle_start']
    genes['gene_circle_end'] = genes['circle_end']
    
    mask = genes['mRNA_name'].isna()
    
    # add gene start and end positions to all expressed canonical gRNAs
    x = gRNAs.merge(genes[~mask][['mO_name', 'cassette_label', 'mRNA_name', 'gene_circle_start', 'gene_circle_end']], on=['mO_name', 'cassette_label', 'mRNA_name'], how='left')
    # for non-expressed canonicals gene position = gRNA position
    x['gene_circle_start'] = x.apply(lambda y: y['circle_start'] if y['expression'] == 'non-expressed' else y['gene_circle_start'], axis=1)
    x['gene_circle_end'] = x.apply(lambda y: y['circle_end'] if y['expression'] == 'non-expressed' or y['gene_circle_end'] is pd.NA else y['gene_circle_end'], axis=1)
    # add gene information on all non-canonical expressed gRNAs 
    x = x.append(genes[mask][['mO_name', 'cassette_label', 'init_seq', 'circle_end', 'circle_start', 'gene_circle_end', 'gene_circle_start', 'expression']])
    # change strand to coding for non-canonicals
    x['strand'] = x['strand'].replace({np.nan:'coding'})
    # change mRNA name to "nc" for non-canonicals
    x['mRNA_name'] = x['mRNA_name'].replace({np.nan:'nc'})
    # create NCBI locus tag for all gRNAs (except non-expressed non-canonicals)
    x['locus_tag'] = x.apply(make_locus_tag, axis=1)
    
    dataframe_out(x, 'tmp.csv')
    
    for mO_name, r in CSB1.items():
        features[mO_name].append((r['start']+1, r['end'], 'misc_feature', {'note':'CSB1'}))
    for mO_name, r in CSB2.items():
        features[mO_name].append((r['start']+1, r['end'], 'misc_feature', {'note':'CSB2'}))
    for mO_name, r in CSB3.items():
        features[mO_name].append((r['start']+1, r['end'], 'misc_feature', {'note':'CSB3'}))
        
    for _, r in cassettes.iterrows():
        features[r['mO_name']].append((r['forward_start']+1, r['reverse_end'], 'misc_structure', {'note':f'gRNA cassette position {r["cassette_label"]}'}))
        features[r['mO_name']].append((r['forward_start']+1, r['forward_end'], 'repeat_region', {'note':'forward', 'rpt_type':'inverted'}))
        features[r['mO_name']].append((r['reverse_start']+1, r['reverse_end'], 'repeat_region', {'note':'reverse', 'rpt_type':'inverted'}))
        
    for _, r in x.iterrows():
        # add the gene (same position as the gRNA by alignment if not expressed)
        if r['strand'] == 'coding':
            features[r['mO_name']].append((r['gene_circle_start']+1, r['gene_circle_end'], 'gene', {'locus_tag':r['locus_tag']}))
        else:
            features[r['mO_name']].append((r['gene_circle_end'], r['gene_circle_start']+1, 'gene', {'locus_tag':r['locus_tag']}))
        # only expressed genes have an initiation sequence
        if r['expression'] == 'expressed':
            if r['strand'] == 'coding':
                features[r['mO_name']].append((r['gene_circle_start']+1, r['gene_circle_start']+init_seq_len, 'misc_feature', {'note':'initiator sequence', 'locus_tag':r['locus_tag']}))
            else:
                features[r['mO_name']].append((r['gene_circle_end'], r['gene_circle_end']-init_seq_len, 'misc_feature', {'note':'initiator sequence', 'locus_tag':r['locus_tag']}))
        # add canonical and non-canonical gRNAs
        if r['mRNA_name'] != 'nc':
            if r['strand'] == 'coding':
                features[r['mO_name']].append((r['circle_start']+1, r['circle_end'], 'ncRNA', {'ncRNA_class':'canonical guide RNA', 'product':r['mRNA_name'], 'note':'guide RNA molecule, specifies the insertion or deletion of nucleotides into the product gene', 'inference':'alignment of minicircle to edited mRNA', 'locus_tag':r['locus_tag']}))
            else:
                features[r['mO_name']].append((r['circle_end'], r['circle_start']+1, 'ncRNA', {'ncRNA_class':'canonical guide RNA', 'product':r['mRNA_name'], 'note':'guide RNA molecule, specifies the insertion or deletion of nucleotides into the product gene', 'inference':'alignment of minicircle to edited mRNA', 'locus_tag':r['locus_tag']}))
        else:
            features[r['mO_name']].append((r['circle_start']+1, r['circle_end'], 'ncRNA', {'ncRNA_class':'non-canonical guide RNA', 'product':'unknown', 'inference':'alignment of RNA-seq to minicircle', 'locus_tag':r['locus_tag']}))
        
    with open(f'{genbank_dir}/minicircles.tbl', 'w') as w:      
        for mO_name, r in features.items():
            w.write(f'>Feature {mO_name}\n')
            # r is a list of feature tuples
            for f in sorted(r, key=itemgetter(0)):
                w.write(f'{f[0]}\t{f[1]}\t{f[2]}\n')
                for k, v in f[3].items():
                    w.write(f'\t\t\t{k}\t{v}\n')
                    
    SeqIO.write(minicircles.values(), f'{genbank_dir}/minicircles.fsa', 'fasta')
    info = \
        '[organism=Trypanosoma brucei brucei]' \
        '[strain=AnTat1.1 90-13]' \
        '[topology=circular]' \
        '[location=mitochondrion]'
        
    cmd = f'tbl2asn -p {genbank_dir} -a s -t minicircles.sbt -V v -j "{info}" -Y comments.txt'
    print(cmd)
    system(cmd)
        
def annotate(minicircles, CSB1, CSB2, CSB3, cassettes, gRNAs, genes, init_seq_len):
    cstrand = {'coding':1, 'template':-1}
    for mO_name, minicircle_record in minicircles.items():
        minicircle_record.description = 'Trypanosoma brucei brucei strain AnTat1.1 90-13'
        minicircle_record.annotations['molecule_type'] = 'DNA'
        # minicircle_record.annotations['accession'] = ''
        # minicircle_record.annotations['version'] = ''
        minicircle_record.annotations['keywords'] = ['WGS']
        minicircle_record.annotations['source'] = 'kinetoplast Trypanosoma brucei brucei'
        minicircle_record.annotations['organism'] = 'Trypanosoma brucei brucei'
        minicircle_record.annotations['topology'] = 'circular'
        minicircle_record.annotations['date'] = str(datetime.datetime.now().strftime("%d-%b-%Y")).upper()
        minicircle_record.annotations['taxonomy'] = ['Eukaryota', 'Euglenozoa', 'Kinetoplastida', 'Trypanosomatidae']
        ref1 = Reference()
        ref1.authors = 'Cooper,S., Wadsworth,E.S., Ochsenreiter,T., Ivens,A., Savill,N.J. and Schnaufer,A.'
        ref1.title = 'Assembly and annotation of the mitochondrial minicircle genome of a differentiation-competent strain of Trypanosoma brucei'
        ref1.journal = 'RNA'
        ref2 = Reference()
        ref2.authors = 'Sinclair Cooper, Elizabeth S. Wadsworth, Achim Schnaufer and Nicholas J. Savill'
        ref2.title = 'Organisation of minicircle cassettes and guide RNA genes of Trypanosoma brucei'
        ref2.journal = 'Submitted'
        minicircle_record.annotations['references'] = [ref1, ref2]

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

        features  = []
        if mO_name in CSB1:
            features += [SeqFeature(FeatureLocation(CSB1[mO_name]['start'], CSB1[mO_name]['end']), type='CSB1')]
        if mO_name in CSB2:
            features += [SeqFeature(FeatureLocation(CSB2[mO_name]['start'], CSB2[mO_name]['end']), type='CSB2')]
        if mO_name in CSB3:
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

    if not config['have transcriptomics']:
        gRNAs['expression'] = 'expressed'

    gRNAs['gene_length'] = gRNAs['gene_rel_end']-gRNAs['gene_rel_start']
    gRNAs['gene_mRNA_end'] = gRNAs['mRNA_end']+gRNAs['rel_pos']
    gRNAs['gene_mRNA_start'] = gRNAs['gene_mRNA_end'] - gRNAs['gene_length']
    gRNAs['sort_pos'] = gRNAs['gene_mRNA_end']

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

        out = []
        j = 1000
        out.append(''.join([str((i//j)%10) if i%j == 0 else ' ' for i in range(1, len(mRNA_seq)+1)]))
        for j in [100, 10, 1]:
            out.append(''.join([str((i//j)%10) if i%(1 if j == 1 else 10) == 0 else ' ' for i in range(1, len(mRNA_seq)+1)]))
        # out.append(''.join(['A' if i == 1 or i == 4 else 'I' if i == 2 else 'U' if i == 3 else '-' for i in mRNA_record['anchor']]))
        # out.append(''.join([str(int(i)) for i in mRNA_record['anchor_count']]))
        out.append(''.join(['M' if (k == 'u' or j != '-') and i == 0 else '-' for i, j, k in zip(mRNA_record['edited'], mRNA_record['deletions'], mRNA_seq)]))
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

                # info = [gRNA['name'], str(gRNA['family_no']), f"{a_type[gRNA['anchor_type']]*gRNA['anchor_len']}"]
                info = []
                # if gRNA['expression'] == 'expressed':
                #     info += ['*']
                info += [gRNA['name'], f"{a_type[gRNA['anchor_type']]*int(gRNA['anchor_len'])}"]
                # info += [gRNA['name']]
                # info += [gRNA['name'], gRNA['family_id']]
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
    # submission(minicircles, CSB1, CSB2, CSB3, cassettes, gRNAs, genes, init_seq_len, genbank_dir)

    if config['output genbank']:
        annotate(minicircles, CSB1, CSB2, CSB3, cassettes, gRNAs, genes, init_seq_len)
        minicircle_list = [v for _, v in sorted(minicircles.items())]
        SeqIO.write(minicircle_list, f'{genbank_dir}/{genbank_file}', 'genbank')

    output_edits(gRNAs, mRNAs, config, alignments_dir)

