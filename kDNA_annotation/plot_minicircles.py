import re
import json
import gzip
import pickle
import pathlib
import datetime
import warnings
import itertools as it
from pprint import pprint
from copy import deepcopy
from collections import Counter, OrderedDict
from operator import itemgetter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation, Reference
# warnings.filterwarnings('ignore')
# plt.style.use('classic')

mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 12
mpl.rcParams['legend.fontsize'] = 12

def get_features(features, no_maxi=True, mRNA=None):
    for mO_name, feature_list in features.items():
        if no_maxi and mO_name == 'Maxicircle':
            continue
        for feature in feature_list:
            if mRNA is None or feature['mRNA_name'] == mRNA:
                yield feature, mO_name

def move_figure(position='top-right'):
    '''
    Move and resize a window to a set of standard positions on the screen.
    Possible positions are:
    top, bottom, left, right, top-left, top-right, bottom-left, bottom-right
    '''

    mgr = plt.get_current_fig_manager()
    px = 1.5*1920
    py = 1.5*1200

    d = 10  # width of the window border in pixels
    if position == 'top':
        # x-top-left-corner, y-top-left-corner, x-width, y-width (in pixels)
        mgr.window.setGeometry(d, 4*d, px - 2*d, py/2 - 4*d)
    elif position == 'bottom':
        mgr.window.setGeometry(d, py/2 + 5*d, px - 2*d, py/2 - 4*d)
    elif position == 'left':
        mgr.window.setGeometry(d, 4*d, px/2 - 2*d, py - 4*d)
    elif position == 'right':
        mgr.window.setGeometry(px/2 + d, 4*d, px/2 - 2*d, py - 4*d)
    elif position == 'top-left':
        mgr.window.setGeometry(d, 4*d, px/2 - 2*d, py/2 - 4*d)
    elif position == 'top-right':
        mgr.window.setGeometry(px/2 + d, 4*d, px/2 - 2*d, py/2 - 4*d)
    elif position == 'bottom-left':
        mgr.window.setGeometry(d, py/2 + 5*d, px/2 - 2*d, py/2 - 4*d)
    elif position == 'bottom-right':
        mgr.window.setGeometry(px/2 + d, py/2 + 5*d, px/2 - 2*d, py/2 - 4*d)

def get_read_depths(filename, minicircles, aligned_gRNAs):
    read_depths = dict([(i, np.zeros(len(m.seq))) for i, m in minicircles.items()])
    with open(filename) as f:
        for line in f:
            x = line.split()
            mO_name = x[0]
            if mO_name[0] == 'n':
                mO_name = 'm'+mO_name[1:]
            if mO_name in read_depths:
                read_depths[mO_name][int(x[1])-1] = int(x[2])

    for gRNA, mO_name in get_features(aligned_gRNAs):
        if mO_name in read_depths:
            gRNA['read_depth'] = max(read_depths[mO_name][gRNA['gRNA_start']:gRNA['gRNA_end']])
    return read_depths

def plot_minicircles(filename, minicircles, cassettes, motifs, aligned_gRNAs, nt_bias_gRNAs, CSB1, CSB2, CSB3, show=True):
    def add_patch(axes, x1, x2, alpha, colour, linewidth=0, fill=True, hatch='', position=None):
        miny, maxy = axes.get_ylim()
        if position == 'top':
            miny = 0.5*(miny+maxy)
        elif position == 'bottom':
            maxy = 0.5*(miny+maxy)
        axes.add_patch(patches.Rectangle((x1, miny), x2-x1, maxy-miny, alpha=alpha, color=colour, linewidth=linewidth, fill=fill, hatch=hatch))

    csb_colour = ('Red', 'DarkRed', 'Purple')
    inter_peak_distance = 40
    ytextpos = 0.25
    fontsize = 12
    scores = {}
    colorcycle = it.cycle(('b', 'g'))
    window = 100
    read_depths = get_read_depths('Small_RNA/read_depth.txt', minicircles, aligned_gRNAs)

    with open(filename) as f:
        for record in json.load(f):
            mO_name, p, offset, length, score = record
            if mO_name not in scores:
                scores[mO_name] = {}
            scores[mO_name][p] = np.array(score)

    for mO_name, m in sorted(minicircles.items()):
        x = np.arange(len(m.seq)-length)
        nc = len(cassettes[mO_name])

        plt.figure(figsize=(9, 6), dpi=figuredpi)
        # plt.figure(figsize=(20, 9.5))
        plt.subplots_adjust(bottom=0.2)
        # move_figure()
        axes = plt.subplot(111)
        axes.set_xlim(0, len(m.seq))
        axes.set_ylim(0.24, 0.42)
        # axes.set_title('minicircle {}'.format(m.id[-3:]))
        axes.set_xlabel('position (nt)')
        axes.set_ylabel('nucleotide bias score')
        # axes.set_yticks([])

        for p in scores[mO_name]:
            axes.plot(x-offset, scores[mO_name][p], lw=0.5)
        # ss = 0.5*(scores[mO_name]['5'][:-inter_peak_distance]+scores[mO_name]['3'][inter_peak_distance:])
        # xs = x[:-inter_peak_distance]
        # axes.plot(xs-offset, ss)

        max_depth = max(read_depths[mO_name])
        # axes.plot(0.25+0.05*read_depths[mO_name]/max_depth, 'k', lw=1.5)
        # xmax = len(read_depths[mO_name])
        # axes.set_xlim(0, xmax)
        # ypos = 0.4
        # cxpos = xmax+20
        position = {'I':0.42, 'II':0.38, 'III':0.34, 'IV':0.30, 'V':0.26 }
        for c in cassettes[mO_name]:
            axes.axvline(c['forward_start'], c='k', lw=1)
            # axes.text(c['forward_start'], ytextpos, c['forward_seq'], size=fontsize)
            axes.axvline(c['reverse_end'],   c='k', lw=1)
            # axes.text(c['reverse_end'], ytextpos-0.005, c['reverse_seq'], ha='right', size=fontsize)
            add_patch(axes, c['forward_start'], c['reverse_end'], 0.2, 'k')
            # ypos = position[c['position']]
            # axes.text(cxpos, ypos, c['position'], color='k', size=fontsize, family='monospace')
            # ypos -= 0.005
            # if len(c['aligned_gRNA_index']) > 0:
            #     gRNA = aligned_gRNAs[mO_name][c['aligned_gRNA_index'][0]]
            #     axes.text(cxpos, ypos, gRNA['mRNA_name'], color='k', size=fontsize, family='monospace')
            #     ypos -= 0.005
            #     axes.text(cxpos, ypos, gRNA['mRNA_seq'], color='k', size=fontsize, family='monospace')
            #     ypos -= 0.005
            #     axes.text(cxpos, ypos, gRNA['pairing'], color='k', size=fontsize, family='monospace')
            #     ypos -= 0.005
            #     axes.text(cxpos, ypos, gRNA['gRNA_seq'], color='k', size=fontsize, family='monospace')
            #     ypos -= 0.005
                # sRNA = gRNA['expressed']
                # if sRNA:
                #     Dstart = gRNA['mRNA_start']-sRNA['mRNA_start']
                #     Dend   = gRNA['mRNA_end']-sRNA['mRNA_end']
                #     if Dstart >= 0:
                #         gstart = 0
                #     else:
                #         gstart = -Dstart
                #     if Dend >= 0:
                #         gend = gRNA['length']-Dend
                #     else:
                #         gend = gRNA['length']
                #     expression = ''.join([' ']*gstart+['-']*(gend-gstart))
                #     axes.text(cxpos, ypos, expression, color='k', size=fontsize, family='monospace')
                #     ypos -= 0.005

            # ypos = ytextpos+0.005
            # for motif_index in c['motif_index']:
            #     ypos += 0.005
            #     motif = motifs[mO_name][motif_index]
            #     add_patch(axes, motif['gRNA_start'], motif['gRNA_end']-1, 0.8, 'c')
                # axes.text(motif['gRNA_start'], ypos, motif['seq'], color='b', size=fontsize)

        for i, csb in enumerate((CSB1, CSB2, CSB3)):
            try:
                add_patch(axes, csb[mO_name]['start'], csb[mO_name]['end']-1, 0.8, csb_colour[i])
            except KeyError:
                pass

        for motif in motifs[mO_name]:
            add_patch(axes, motif['gRNA_start'], motif['gRNA_end']-1, 0.8, 'b')

        ypos = 0.41
        for gRNA in sorted(aligned_gRNAs[mO_name], key=itemgetter('gRNA_start')):
            if gRNA['cassette_pos'] == 'Orphan':
                add_patch(axes, gRNA['gRNA_start'], gRNA['gRNA_end'], 0.2, 'y')
            else:
                add_patch(axes, gRNA['gRNA_start'], gRNA['gRNA_end'], 0.2, 'r')
                axes.text(gRNA['gRNA_start'], ypos, gRNA['name'], size=fontsize)
                ypos -= 0.01

        for gRNA in nt_bias_gRNAs[mO_name]:
            axes.axvline(gRNA['gRNA_start'], c='r', lw=1)
            axes.axvline(gRNA['gRNA_end'],   c='r', lw=1)
        if show:
            plt.show()

##########################################################################################
figuredpi = 300
gaps = {True:'gaps', False:'nogaps'}
expds = {True:'all', False:'exp'}
analysis = {'allow_gaps':False, 'allow_non_expressed':True}

date = datetime.datetime.now().strftime("%Y-%m-%d")
load_date = '2018-12-13'

ext = '{}_{}'.format(gaps[analysis['allow_gaps']], expds[analysis['allow_non_expressed']])

# load in features from first pass
with gzip.open('Features/all_pass_2_pickled_{}_{}.txt.gz'.format(ext, load_date), 'rb') as f:
    minicircles = pickle.load(f)
    mRNAs = pickle.load(f)
    CSB1 = pickle.load(f)
    CSB2 = pickle.load(f)
    CSB3 = pickle.load(f)
    cassettes = pickle.load(f)
    aligned_gRNAs = pickle.load(f)
    nt_bias_gRNAs = pickle.load(f)
    small_RNAs = pickle.load(f)
    motifs = pickle.load(f)
    editing_groups = pickle.load(f)

# Plotting
# plot_minicircles('json/aligned_gRNAs_scores.json', minicircles, cassettes, motifs, aligned_gRNAs, nt_bias_gRNAs)

for mO_name in ['mO_170']:
    plot_minicircles('json/aligned_gRNAs_scores.json', {mO_name:minicircles[mO_name]}, cassettes, motifs, aligned_gRNAs, nt_bias_gRNAs, CSB1, CSB2, CSB3, show=False)
    # plt.savefig('Paper/Figures/{}_{}.pdf'.format(mO_name, date))
    plt.show()

