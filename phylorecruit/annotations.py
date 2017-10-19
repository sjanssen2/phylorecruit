import sys
import pandas as pd
import urllib.request
import os
from ggmap.snippets import cache


def create_mapping_uniprot_90to50(file_uniprotmapping, err=sys.stderr):
    # this function takes some time!
    # without decompression:        Wall time: 5min 32s
    # with decompression in pandas: Wall time: 8min 47s

    # make sure to use the right path
    file_uniprotmapping = os.path.abspath(file_uniprotmapping)

    # if huge file does not exist, attemp to download it from uniprot
    if not os.path.exists(file_uniprotmapping):
        if not os.path.exists(file_uniprotmapping+'.gz'):
            urllib.request.urlretrieve(
                ('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release'
                 '/knowledgebase/idmapping/idmapping_selected.tab.gz'),
                file_uniprotmapping+'.gz')

    if not os.path.exists(file_uniprotmapping+'.map'):
        # actual reading of the file
        if err is not None:
            sys.stderr.write('1) parsing file.\n')
        x = pd.read_csv(file_uniprotmapping+'.gz',
                        sep='\t',
                        header=None,
                        usecols=[8, 9],
                        names=['uniprot90', 'uniprot50'],
                        dtype=str,
                        compression='gzip',
                        engine='c').dropna()

        # remove the 'uniprotXX_' prefix
        if err is not None:
            sys.stderr.write('2) removing prefix.\n')
        for c in x.columns:
            x[c] = x[c].apply(lambda y: y[9:])

        # filter to only those instances where the cluster name between
        # uniprot50 is actually different from uniprot90
        if err is not None:
            sys.stderr.write('3) filter identities.\n')
        x = x[x['uniprot50'] != x['uniprot90']]

        # drop duplicate information
        if err is not None:
            sys.stderr.write('4) drop duplicates.\n')
        x.drop_duplicates(inplace=True)

        if err is not None:
            sys.stderr.write('5) write output file.\n')
        x.to_csv(file_uniprotmapping+'.map', sep='\t', index=False)
    else:
        if err is not None:
            sys.stderr.write('1) reading cached mapping.\n')
        x = pd.read_csv(file_uniprotmapping+'.map', sep='\t', dtype=str)
    x = x.set_index('uniprot90')
    return x


@cache
def obtain_counts_ecnumbers(genomes, err=sys.stderr):
    counts = dict()
    coordinates = dict()
    if err is not None:
        err.write('parsing EC numbers for %i genomes: ' % genomes.shape[0])
    for genome, row in genomes.iterrows():
        if err is not None:
            err.write('.')
        counts[genome] = dict()
        coordinates[genome] = dict()
        f = open(row['path'], 'r')
        for line in f:
            if 'EC_number=' in line:
                fields = line.split('\t')
                ecs = fields[8].split(';')[0].split('=')[1].split(',')
                for ec in ecs:
                    if ec not in counts[genome]:
                        counts[genome][ec] = 0
                    counts[genome][ec] += 1
                    if ec not in coordinates[genome]:
                        coordinates[genome][ec] = []
                    coordinates[genome][ec].append((fields[0],
                                                    fields[3],
                                                    fields[4]))
        f.close()

    if err is not None:
        err.write(' done.\n')

    # collapse EC numbers at all four levels
    results = dict()
    counts = pd.DataFrame(counts).fillna(0)
    for level in list(range(1, 5)):
        counts['EC'] = list(map(lambda x: '.'.join(x.split('.')[:level]),
                                counts.index))
        results['level_%i' % level] = counts.groupby('EC').sum()

    return results, coordinates


@cache
def obtain_counts_uniprot(genomes, mapping, err=sys.stderr):
    counts = dict()
    coordinates = dict()
    if err is not None:
        err.write('parsing uniprot IDs for %i genomes: ' % genomes.shape[0])
    for genome, row in genomes.iterrows():
        if err is not None:
            err.write('.')
        counts[genome] = dict()
        coordinates[genome] = dict()
        file_diamond = '/'.join(row['path'].split('/')[:-1])+'/tmp/diamond.hit'
        if os.path.exists(file_diamond):
            hits = pd.read_csv(file_diamond, sep='\t', dtype=str, header=None,
                               usecols=[0, 2], names=['contig', 'uniprot'])
            hits['identity'] = hits['uniprot'].apply(lambda x: x[6:8])
            hits['uniprot'] = list(map(lambda x: x[9:], hits['uniprot']))

            for _, row in hits.iterrows():
                _id = None
                # uniprot hits against 50% clusters don't need to be mapped
                # at all
                if row['identity'] == '50':
                    _id = row['uniprot']
                else:
                    try:
                        # map 90% cluster names to 50% cluster names if they
                        # are not the same, i.e. are in the mapping table
                        _id = mapping.loc[row['uniprot'], 'uniprot50']
                    except KeyError:
                        # otherwise name of 50% and 90% cluster are identical
                        _id = row['uniprot']
                if _id not in counts[genome]:
                    counts[genome][_id] = 0
                counts[genome][_id] += 1
                if _id not in coordinates[genome]:
                    coordinates[genome][_id] = []
                coordinates[genome][_id].append((row['contig'], 0, 0))

    if err is not None:
        err.write(' done.\n')
    return {'uniprot50': pd.DataFrame(counts).fillna(0)}, coordinates
