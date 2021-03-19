import re
import gzip
import yaml
import pickle
import numpy as np
import pandas as pd
from Bio import SeqIO
from copy import copy
from Bio.Seq import Seq

pd.options.display.max_colwidth = 100

def complement(seq):
    conv = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', '-':'-'}
    try:
        return ''.join([conv[i] for i in seq])
    except:
        # non-standard base found
        return None

def pairs(mRNA, gRNA):
    pairings = {'GC':'|', 'CG':'|', 'AU':'|', 'UA':'|', 'GU':':', 'UG':':'}
    return ''.join([pairings[mb+gb] if mb+gb in pairings else '-' if '-' in mb+gb else '.' for mb, gb in zip(mRNA, gRNA)])

def get_mRNAs(insertion_file, deletion_file):
    #edits:   ACGTt
    #DNA_seq: ACGT
    #seq:     ACGUu

    mRNA_name = ''
    deletions = {}
    with open(deletion_file) as f:
        for line in f:
            if line.startswith('>'):
                mRNA_name = line[1:-1]
            else:
                deletions[mRNA_name] = line[:-1]

    insertions = {}
    with open(insertion_file) as f:
        for line in f:
            if line.startswith('>'):
                mRNA_name = line[1:-1]
            else:
                insertions[mRNA_name] = line[:-1]

    def add_N(seq):
        # correct length of truncated sequence to be able to translate it
        l = len(seq) % 3
        if l == 0:
            return seq
        elif l == 1:
            return seq+'NN'
        else:
            return seq+'N'

    mRNAs = {}
    for name, edits in sorted(insertions.items()):
        #determine open reading frame as the one that gives the longest protein sequence between 2 stop codons
        mRNA_seq = edits.replace('t', 'U').replace('T', 'U')
        max_lengths = []
        for orf in range(3):
            max_lengths.append(max( [len(j) for j in re.split(r'\*', str(Seq(add_N(mRNA_seq[orf:])).translate(table=4)))] ))
        orf = max_lengths.index(max(max_lengths))

        mRNAs[name] = {'orf'      :orf,
                       'seq'      :edits.replace('t', 'u').replace('T', 'U'),
                       'edits'    :edits,
                       'length'   :len(edits),
                       'deletions':deletions[name],
                       'DNA_seq'  :edits.replace('t', 'T'),
                       'translate':Seq(add_N(mRNA_seq[orf:])).translate(table=4)
                       }
    print(f'Number of mRNAs: {len(mRNAs)}')
    return mRNAs

def get_minicircles(filename, file_format='fasta'):
    m = SeqIO.to_dict(SeqIO.parse(filename, file_format))
    print(f'Number of minicircles: {len(m)}')
    return m

def get_maxicircle(filename, file_format='fasta'):
    m = SeqIO.to_dict(SeqIO.parse(filename, file_format))
    return m

def pickle_save(x, file):
    with open(file, 'wb') as f:
        pickle.dump(x, f)
    print(f'file {file} written')

def pickle_load(file):
    with open(file, 'rb') as f:
        print(f'file {file} loaded')
        return pickle.load(f)

def gzip_pickle_save(x, file):
    with gzip.open(file, 'wb') as f:
        pickle.dump(x, f)
    print(f'file {file} written')

def gzip_pickle_load(file):
    with gzip.open(file, 'rb') as f:
        print(f'file {file} loaded')
        return pickle.load(f)

def dataframe_out(x, file, index=True):
    with open(file, 'w') as f:
        f.write(x.to_string(index=index))
    print(f'file {file} written')

def get_sequence(gRNA, up, down, minicircles):
    if gRNA['strand'] == 'coding':
        s = max(0, gRNA['circle_start']-up)
        e = gRNA['circle_start']+down
    else:
        s = max(0, gRNA['circle_end']-down)
        e = gRNA['circle_end']+up
    return str(minicircles[gRNA['mO_name']].seq)[s:e]

def load_config(file):
    with open(file) as f:
        print(f'file {file} loaded')
        return yaml.safe_load(f)

def get_directories(config):
    project = config['project']
    in_dir = f"{project}/{config['in directory']}" 
    work_dir = f"{project}/{config['working directory']}" 
    annotation_dir = f"{project}/{config['annotation directory']}" 
    meme_dir = f"{project}/{config['meme directory']}" 
    return in_dir, work_dir, annotation_dir, meme_dir

def identify_anchors(gRNAs, mRNAs, filter):
    """
    For each gRNA we identify its maximum anchor as the longest, 5’-most 
    stretch of Watson-Crick basepairs. An anchor is considered an extender 
    if there are any insertions or deletions in its first 6 basepairs, an 
    initiator if there are no insertions or deletions in its first 6 basepairs.
    """

    # dictionary for creating an anchor dataframe to join onto gRNAs dataframe
    columns = ('anchor_type',)
    anchor = dict([(c, []) for c in columns])
    index = []

    anchor_regex = re.compile(r'\|{{{},}}'.format(filter['min_anchor_length']))
    guiding_regex = re.compile(r'[\|:\.]+')

    for mRNA_name, mRNA in mRNAs.items():
        # track prior editing and anchors on the mRNA
        mRNA['anchor_count'] = np.zeros(len(mRNA['seq']))
        mRNA['anchor_len'] = np.zeros(len(mRNA['seq']))
        mRNA['edited']  = np.zeros(len(mRNA['seq']))

        # all edits (insertions and deletions along) an mRNA
        edits = [True if i == 'u' or j != '-' else False for i, j in zip(mRNA['seq'], mRNA['deletions'])]

        # get all gRNAs for this mRNA
        g = copy(gRNAs[gRNAs['mRNA_name'] == mRNA_name])
        # and sort on mRNA_end and length
        g = g.sort_values(['mRNA_end', 'length'], ascending=[False, True])

        for idx, gRNA in g.iterrows():
            # find position of minimal anchor on gRNA and the mRNA sequence of this anchor
            match = anchor_regex.search(gRNA['pairing'][::-1])
            if match is None:
                gRNA['anchor_len'] = None
                continue
            # find position of anchor in alignment pairing
            a_start = match.start(0)
            a_end   = match.end(0)
            # maximal anchor on mRNA
            max_a_slice = slice(gRNA['mRNA_end']-a_end, gRNA['mRNA_end']-a_start)
            # minimal anchor on mRNA (ie min_anchor_length)
            min_a_slice = slice(gRNA['mRNA_end']-a_start-filter['min_anchor_length'], gRNA['mRNA_end']-a_start)
            # get mRNA sequence of min and max anchors
            min_seq  = edits[min_a_slice]
            max_seq  = edits[max_a_slice]
            # min_seq  = mRNA['seq'][min_a_slice]
            # max_seq  = mRNA['seq'][max_a_slice]
            min_prior_edits = mRNA['edited'][min_a_slice]
            max_prior_edits = mRNA['edited'][max_a_slice]

            # determine the type of anchor of this gRNA
            # if 'u' not in min_seq:
            if True not in min_seq:
                if np.sum(min_prior_edits) == 0:
                    a_type = 'initiator'
                else:
                    a_type = 'extenderB' # based on min_anchor_length could be an initiator or an extender
            else:
                for m, e in zip(min_seq, min_prior_edits):
                    # if m == 'u' and e == 0:
                    if m and e == 0:
                        a_type = 'unanchored'
                        break
                else:
                    a_type = 'extenderA'

            # update position of end of anchor based on anchor type and position of prior edits
            if a_type == 'extenderA':
                # Normal anchor using minimum WC base-pairing
                # ie anchor is created by prior insertions
                # pos = np.where(mRNA['edited'][max_a_slice] > 0)[0][0]
                # find 5'-most position of last insertion of prior gRNA edits 
                for p, (m, e) in enumerate(zip(max_seq[::-1], max_prior_edits[::-1])):
                    # if e == 0 and m == 'u':
                    if e == 0 and m:
                        pos = len(max_seq)-p
                        a_end -= pos
                        break
                a_value = 1
            elif a_type == 'extenderB':
                # Normal anchor using minimum WC base-pairing
                # ie anchor is created by prior insertions
                # pos = np.where(mRNA['edited'][max_a_slice] > 0)[0][0]
                # find 5'-most position of last insertion of prior gRNA edits 
                for p, (m, e) in enumerate(zip(max_seq[::-1], max_prior_edits[::-1])):
                    # if e == 0 and m == 'u':
                    if e == 0 and m:
                        pos = len(max_seq)-p
                        a_end -= pos
                        break
                a_value = 4
            elif a_type == 'initiator':
                # # Initiator anchor of at least minimum length
                # # trim maximum anchor to first 'u' if it exists
                # match = re.search(r'\d', mRNA['deletions'][max_a_slice][::-1])
                # # search for deletions in anchor (this only seems to happen in A6_v1)
                # if match:
                #     dpos = match.start(0)
                # else:
                #     dpos = len(max_seq)
                # # find first insertion in anchor
                # # epos = max_seq[::-1].find('u')
                # epos = max_seq[::-1].find('u')
                # if epos == -1:
                #     epos = len(max_seq)
                # # take the minimum position of these two
                # pos = min(dpos, epos)
                try:
                    pos = max_seq.index(True)
                except ValueError:
                    pos = len(max_seq)
                a_end -= len(max_seq)-pos
                a_value = 2
            else:
                # otherwise gRNA has no edits to anchor it, unattached anchor
                # in analysis.py this is considered an extender
                # position end of anchor at minimum anchor length away form start of anchor
                a_end = a_start+filter['min_anchor_length']
                a_value = 3

            # add anchor region to mRNA and gRNA
            a_slice = slice(gRNA['mRNA_end']-a_end, gRNA['mRNA_end']-a_start)
            # a_slice = slice(gRNA['mRNA_end']-a_start-filter['min_anchor_length'], gRNA['mRNA_end']-a_start)
            mRNA['anchor_count'][a_slice] += 1
            x = mRNA['anchor_len'][a_slice]
            np.place(x, (x==0) | (x==3), a_value)

            index.append(idx)
            anchor['anchor_type'].append(a_type)

            # add edited region to mRNA
            match = guiding_regex.match(gRNA['pairing'][::-1][a_end:])
            if match:
                e_start = match.start(0)+a_end
                e_end   = match.end(0)+a_end
                mRNA['edited'][gRNA['mRNA_end']-e_end:gRNA['mRNA_end']-e_start] += 1
        mRNA['anchor_count'] = np.clip(mRNA['anchor_count'], 0, 9)
    return gRNAs.join(pd.DataFrame(anchor, index=index)), mRNAs

def orphan_position(gRNA, cassettes):
    mO_name = gRNA['mO_name']
    # extract cassettes on the minicircle encoding this gRNA
    cassettes = cassettes[cassettes['mO_name'] == mO_name]
    # No cassettes on this gRNA (eg maxicircle)
    if len(cassettes) == 0:
        return None

    # order cassettes along minicircle
    ordered = cassettes.sort_values('forward_start')
    # lists of start of forward repeat and ends of reverse repeat
    cassette_starts = ordered['forward_start'].values
    cassette_ends = ordered['reverse_end'].values

    if gRNA['circle_start'] < cassette_starts[0]:
        c = ordered.iloc[0]
        return f'before {c["cassette_label"]}'

    if gRNA['circle_start'] > cassette_ends[-1]:
        c = ordered.iloc[-1]
        return f'after {c["cassette_label"]}'

    # look for gRNAs the end of one cassette and the start of the next
    for idx, (x1, x2) in enumerate(zip(cassette_ends[:-1], cassette_starts[1:])):
        if x1 < gRNA['circle_start'] < x2:
            c1 = ordered.iloc[idx]
            c2 = ordered.iloc[idx+1]
            return f'between {c1["cassette_label"]} and {c2["cassette_label"]}'
    
    # otherwise gRNA is within a cassette
    return None

def trim_mismatches(gRNA):
    """ 
        trim off 3' mismatches
        1. first trim to prevent long alignments past end of normal expressed gRNAs
        2. second trim to mismatches close to end of gRNA
    """
    pairing = gRNA['pairing']

    # 1. index of right-most mismatch before index -40
    # if no MM is found trim will equal 0 
    trim = pairing.rfind('.', 0, -40)+1

    # 2. 
    while '.' in pairing[trim:]:
        mm = pairing.find('.', trim)
        # if distance to next mismatch is less than 4 then trim
        # otherwise stop trimming
        if mm-trim < 4:
            # trim += mm+1
            trim = mm+1
        else:
            break

    # trim alignment information
    if trim > 0:
        if gRNA['strand'] == 'coding':
            gRNA['circle_end'] -= trim
        else:
            gRNA['circle_start'] += trim
        gRNA['mRNA_start'] += trim
        gRNA['length']     -= trim
        gRNA['mRNA_seq']    = gRNA['mRNA_seq'][trim:]
        gRNA['gRNA_seq']    = gRNA['gRNA_seq'][trim:]
        gRNA['pairing']     = gRNA['pairing'][trim:]
        gRNA['mismatches']  = gRNA['pairing'].count('.')
    return gRNA

def collapse(gRNAs, minicircles, mRNAs, post_cassette_label=True):
    if len(gRNAs) == 1:
        return gRNAs
    else:
        strand = gRNAs['strand'].values[0]
        circle_name = gRNAs['mO_name'].values[0]
        mRNA_name = gRNAs['mRNA_name'].values[0]
        # merge gRNA sequences
        min_circle_start = gRNAs['circle_start'].min()
        max_circle_end   = gRNAs['circle_end'].max()
        if post_cassette_label:
            min_rel_start = gRNAs['rel_start'].min()
        else:
             min_rel_start = None
        if strand == 'coding':
            gRNA_seq = str(minicircles[circle_name].seq[min_circle_start:max_circle_end]).replace('T', 'U')[::-1]
        else:
            gRNA_seq = complement(str(minicircles[circle_name].seq[min_circle_start:max_circle_end])).replace('T', 'U')
        # merge mRNA sequences
        min_mRNA_start = gRNAs['mRNA_start'].min()
        max_mRNA_end   = gRNAs['mRNA_end'].max()
        mRNA_seq = mRNAs[mRNA_name]['edits'][min_mRNA_start:max_mRNA_end].replace('t', 'u').replace('T', 'U')

        # indicative of multiple gRNAs editing same mRNA but at different positions so keep all
        # to prevent creating very long gRNAs
        if len(gRNA_seq) != len(mRNA_seq):
            return gRNAs

        # construct merged alignment of collapsed alignments
        gRNA = {}
        gRNA['mO_name']      = circle_name
        if post_cassette_label:
            gRNA['cassette_label'] = gRNAs['cassette_label'].values[0]
        gRNA['strand']       = strand
        gRNA['length']       = max_circle_end-min_circle_start
        if post_cassette_label:
            gRNA['rel_start']    = min_rel_start
        gRNA['circle_start'] = min_circle_start
        gRNA['circle_end']   = max_circle_end
        gRNA['mRNA_name']    = mRNA_name
        gRNA['product']      = gRNAs['product'].values[0]
        gRNA['mRNA_start']   = min_mRNA_start
        gRNA['mRNA_end']     = max_mRNA_end
        gRNA['mRNA_seq']     = mRNA_seq
        gRNA['gRNA_seq']     = gRNA_seq
        gRNA['pairing']      = pairs(mRNA_seq.upper(), gRNA_seq)
        gRNA['mismatches']   = gRNA['pairing'].count('.')
        gRNA['anchor_len']   = max(gRNAs['anchor_len'])
        return pd.DataFrame(gRNA, index=[0])

def drop_smaller(gRNAs):
    if len(gRNAs) == 1:
        return gRNAs
    else:
        return gRNAs.loc[[gRNAs['length'].astype(int).idxmax()]]

