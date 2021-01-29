import re
import sys
import datetime
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter as cC

mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['legend.fontsize'] = 14

bg_colour = 'LightGray'
all_colours = {'CSB1':'Red', 'CSB2':'DarkRed', 'CSB3':'Purple', 'CASSETTE':'Blue', 'small_RNA':'Yellow'}

def get_minicircles(minicircles, filename):
    gbk_feature_list = all_colours.keys()
    posregex = re.compile(r'(\d+)\.{2}(\d+)')
    seqregex = re.compile('[acgt]*')
    maxlength = 0

    with open(filename) as o:
        for line in o:
            if line.startswith('LOCUS'):
                name = line.split()[1]
                x = []
            for j, f in enumerate(gbk_feature_list):
                if f in line:
                    pos = posregex.search(line)
                    try:
                        x.append([int(pos.group(1)), int(pos.group(2)), f])
                    except AttributeError:
                        pass
            if line.startswith('ORIGIN'):
                s = []
                while line[0] != '/':
                    s += seqregex.findall(line)
                    line = next(o)
                sequence = ''.join(s)
                n = len(sequence)
                maxlength = max(maxlength, n)
                minicircles[name] = {'sequence':sequence, 'features':x, 'length':n}
    return maxlength

def get_cassette_locations(mO_filename):
    mOs = SeqIO.to_dict(SeqIO.parse(mO_filename, 'genbank'))

    positions = {}
    for mO, record in mOs.items():
        for f in record.features:
            if f.type == 'CASSETTE':
                pos = f.qualifiers['position'][0]
                if pos not in positions:
                    positions[pos] = []
                positions[pos].append(f.location.start+(f.location.end-f.location.start)/2)

    for pos in positions:
        print(pos, np.mean(positions[pos]))
        print(pos, np.std(positions[pos]))
    return positions

def feature_index_position(feature_list, feature):
    for i, f in enumerate(feature_list):
        if feature in f[2]: return i, f[0], f[1]
    return None, None, None
    
def plot_features(minicircles, maxlength, locations, align='CSB1', position='start', colours=None, centre=True):
    sort_by_length = sorted(minicircles.items(), key=lambda x: x[1]['length'])

    #find start of aligned feature
    pos, l1, l2 = [], [], []
    for minicircle in sort_by_length:
        structure = minicircles[minicircle[0]]
        index, p0, p1 = feature_index_position(structure['features'], align)
        if position == 'end':
            pos.append(p1)
        else:
            pos.append(p0)
        if index is None:
            l1.append(None)
            l2.append(None)
        else:
            l1.append(pos[-1]-1)
            l2.append(structure['length']-pos[-1])

    row = 0
    L1 = max(l1)
    img = np.ones((len(pos)-pos.count(None), maxlength, 3))*cC.to_rgb('White')

    for i, minicircle in enumerate(sort_by_length):
        if pos[i] is None: 
            continue
        name = minicircle[0]
        structure = minicircles[name]
        features = structure['features']
        l = structure['length']
        d = L1-pos[i]+1
        d2 = l/2
        mid = maxlength/2
        img[row, d:l+d] = cC.to_rgb(bg_colour)
        for f in features:
            img[row, (f[0]+d-1):(f[1]+d)] = cC.to_rgb(colours[f[2]])
        x = np.array(list(structure['sequence']))
        i = np.where(x=='a')[0]
        i = i[(i>l-70) & (i<l-20)]
        img[row, np.add(d, i)] = cC.to_rgb('DarkGrey')
        row += 1

    fig = plt.imshow(img, interpolation='none')
    plt.xlabel('Position (nt)')
    plt.ylabel('Minicircle')
    plt.yticks([0, 99, 199, 299, 397], [1, 100, 200, 300, 398])
    for l in locations:
        plt.annotate(l, (np.median(locations[l]), 0), size=14)
    plt.savefig(f'{new_paper}/structure_{date}.png', dpi=figuredpi, bbox_inches='tight')
    plt.show()

project = 'Antat1.1'
load_date = '2020-12-28'

directory = f'{project}/Annotation/{load_date}/Minicircles'
mO_filename = f'{directory}/annotated_mOs_{load_date}.gbk'
new_paper = f'{project}/Paper/Figures'
figuredpi = 300
date = datetime.datetime.now().strftime("%Y-%m-%d")
minicircles = {}
maxlength = get_minicircles(minicircles, mO_filename)
cassette_locations = get_cassette_locations(mO_filename)
plot_features(minicircles, maxlength, cassette_locations, align='CSB1', position='start', colours=all_colours, centre=False)

