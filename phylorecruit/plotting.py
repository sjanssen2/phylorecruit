from ete3 import Tree, SeqMotifFace
import pandas as pd
import seaborn as sns


availColors = sns.color_palette('Paired', 12) + sns.color_palette('Dark2', 12) + sns.color_palette('Pastel1', 12)


def boxes_presabs(tree, cluster, coords, gene_size=60):
    """
    Returns
    -------
    (boxes, size), where
    - boxes is a dict (key=genomes) holding lists of annotation boxes for ete
    - and size is the maximal width of all annotations.
    """
    boxes = dict()
    colormap = dict()
    annot_position = dict()
    for node in tree.traverse():
        if node.is_leaf():
            genome = node.name
            hits = dict()
            for ec in cluster:
                if ec in coords[genome]:
                    for contig, start, end in coords[genome][ec]:
                        if ec not in annot_position:
                            annot_position[ec] = len(annot_position)
                        pos = annot_position[ec]
                        if ec not in hits:
                            hits[ec] = 0
                        hits[ec] += 1
                        #hits.append({'annotation': ec})

            for ec in hits.keys():
                pos = annot_position[ec]
                if ec not in colormap:
                    colormap[ec] = '#%02x%02x%02x' % tuple(map(lambda x: int(255*x), availColors[len(colormap)]))
                color = colormap[ec]
                if genome not in boxes:
                    boxes[genome] = []
                boxes[genome].append([gene_size*pos, gene_size*(1+pos), "[]", None, 10, 'black', color, "arial|8|black|%s: %i" % (ec, hits[ec])])

            if genome in boxes:
                boxes[genome] = sorted(boxes[genome], key=lambda x: x[0])
    return boxes, len(annot_position) * gene_size


def boxes_position(tree, cluster, coords, genomes, size=1000, min_gene_size=5):
    """
    Returns
    -------
    (boxes, size), where
    - boxes is a dict (key=genomes) holding lists of annotation boxes for ete
    - and size is the maximal width of all annotations.
    """
    boxes = dict()
    colormap = dict()
    for node in tree.traverse():
        if node.is_leaf():
            genome = node.name
            hits = []
            for ec in cluster:
                if ec in coords[genome]:
                    for contig, start, end in coords[genome][ec]:
                        hits.append({'annotation': ec,
                                     'contig': contig,
                                     'start': int(start),
                                     'end': int(end)})
            if len(hits) > 0:
                boxes[genome] = []
                hits = pd.DataFrame(hits)
                offset = 0
                lens = list(map(int, str(genomes.loc[genome, 'contig_lengths']).split(';')))
                names = genomes.loc[genome, 'contig_names'].split(';')
                genome_size = sum(lens)
                for contig, contig_length in sorted(zip(names, lens)):
                    idx = hits[hits['contig'] == contig].index

                    # correct chromosomal coordinates to genome wide coordinates
                    hits.loc[idx, 'start'] += offset
                    hits.loc[idx, 'end'] += offset

                    # map genomic coordinates into given drawing space
                    hits.loc[idx, 'start'] = ((size / genome_size) * hits.loc[idx, 'start']).astype(int)
                    hits.loc[idx, 'end'] = ((size / genome_size) * hits.loc[idx, 'end']).astype(int)

                    # ensure that genes have a minimal length to be visible in the plot
                    hits.loc[idx, 'end'] = hits.loc[idx, :].apply(lambda row: max(row['start']+min_gene_size, row['end']), axis=1)

                    # ete3 boxes for annotations
                    for i, annot in hits.loc[idx,:].iterrows():
                        if annot['annotation'] not in colormap:
                            colormap[annot['annotation']] = '#%02x%02x%02x' % tuple(map(lambda x: int(255*x), availColors[len(colormap)]))
                        color = colormap[annot['annotation']]
                        boxes[genome].append([hits.loc[i, 'start'], hits.loc[i, 'end'], "[]", None, 10, 'black', color, "arial|8|black|%s" % annot['annotation']])

                    offset += contig_length
    return (boxes, size)


def plot_tree(file_tree, cluster, coords, genomes, size=1000, file_result=None, kind='presAbs'):
    tree = Tree(file_tree, format=1)
    if kind == 'presAbs':
        boxes, size = boxes_presabs(tree, cluster, coords)
    else:
        boxes, size = boxes_position(tree, cluster, coords, genomes, size=size)

    for node in tree.traverse():
        if node.is_leaf():
            rightmost = 0
            allboxes = []
            if (node.name in boxes) and (len(boxes[node.name]) > 0):
                allboxes.extend(boxes[node.name])
                rightmost = boxes[node.name][-1][1]
            if rightmost < size:
                color = 'lightgray'
                if rightmost > 0:
                    color = 'black'
                allboxes.append([rightmost, size, "line", None, 10, color, None, "arial|5|black|"])
            seqFace = SeqMotifFace(seq=None, motifs=allboxes)
            node.add_face(seqFace, 0, "aligned")

    if file_result is None:
        return tree.render("%%inline")
    else:
        tree.render(file_result)
