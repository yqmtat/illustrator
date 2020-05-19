import csv, codecs
from Bio import Entrez
from Bio import SeqIO
import os
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram

Entrez.email = "example@example.com"
input = 'data.csv'
output = 'phage.pdf'
CLASS1=['NC_005856.1','NC_000866.4','NC_005859.1']
CLASS2=['LT960609.1','LT960607.1','NC_001604.1']

def readcsv(file):
    with open(file, 'r') as f:
        reader = csv.reader(f)
        data = list(reader)
        new_data = []
        for _ in data[1:]:
            new_data.append(
                [_[0], _[1], _[2].split('-')[0], _[2].split('-')[1], _[3].split('-')[0], _[3].split('-')[1]])
    return new_data


new_data = readcsv(input)

GenomeId = []
for i in range(len(new_data)):
    contig = []
    for pos in [1, 2]:
        Id = new_data[i][pos - 1].split(' ')[0]
        num_start = int(new_data[i][2 * pos])
        num_end = int(new_data[i][2 * pos + 1])
        if Id not in GenomeId:
            try:
                seq = SeqIO.parse(Id + ".fasta", "fasta")
            except FileNotFoundError:
                hd = Entrez.efetch(db="nucleotide", id=Id, rettype='fasta')
                print('下载' + Id + '.fasta')
                seq = SeqIO.read(hd, 'fasta')
                fw = open(Id + '.fasta', 'w')
                SeqIO.write(seq, fw, 'fasta')
                fw.close()
                os.getcwd()
            GenomeId.append(Id)
        for seq_record in SeqIO.parse(Id + ".fasta", "fasta"):
            print("导入" + Id + '.fasta')
            contig.append(seq_record.seq[num_start:num_end])

    if contig[0] == contig[1]:
        new_data[i].append(contig[0])
        new_data[i].append(1)
    elif contig[0] == contig[1].reverse_complement():
        new_data[i].append(contig[0])
        new_data[i].append(-1)
    else:
        new_data[i].append(contig[0])
        new_data[i].append(contig[1])


def data_write_csv(file_name, datas):
    file_csv = codecs.open(file_name, 'w+', 'utf-8')
    writer = csv.writer(file_csv)
    for data in datas:
        writer.writerow(data)
    print('共有序列下载完成！')


def GenomeMap(file, GenomeId, grid=10000, cross=True):
    # print(GenomeId)
    gd_diagram = GenomeDiagram.Diagram('phages')
    with open(file, 'r') as f:
        reader = csv.reader(f)
        data = list(reader)

    records = []
    ref = {}
    for Id in GenomeId:
        try:
            record = SeqIO.read(Id + ".gb", "genbank")
        except FileNotFoundError or IOError or ValueError:
            hd = Entrez.efetch(db="nucleotide", id=Id, rettype='gb', retmode="text")
            record = SeqIO.read(hd, 'genbank')
            fw = open(Id + '.gb', 'w')
            SeqIO.write(record, fw, 'genbank')
            fw.close()
            os.getcwd()
        for i in SeqIO.parse(Id + ".gb", "genbank"):
            ref[Id] = i.annotations['keywords']
        records.append(record)

    feature_sets = {}
    max_len = 0
    for i, record in enumerate(records):
        max_len = max(max_len, len(record))
        gd_track_for_features = gd_diagram.new_track(
            5 - 2 * i,
            name=record.description,
            greytrack=True,
            greytrack_fontsize=16,
            greytrack_labels=1,
            largetick=True,
            smalltick=True,
            scale_ticks=True,
            scale_largeticks=0.5,
            scale_smallticks=0.1,
            scale_largetick_interval=grid,
            scale_smalltick_interval=grid / 20,
            scale_largetick_labels=True,

            start=0,
            end=len(record),
        )
        assert record.name not in feature_sets
        feature_sets[record.id] = gd_track_for_features.new_set()

    for crosslink in data:
        if not cross:
            break
        set_X = feature_sets[crosslink[0].split(' ')[0]]
        set_Y = feature_sets[crosslink[1].split(' ')[0]]
        # 手动划分连接类型时使用
        # score = 100
        # try:
        #     if crosslink[7] == 1 or crosslink[7] == -1:
        #         score = 100
        # except TypeError:
        #     score = 50
        if crosslink[0].split(' ')[0] in CLASS1 and crosslink[1].split(' ')[0] in CLASS1:
            color = colors.linearlyInterpolatedColor(colors.green, colors.yellow, 0, len(GenomeId),
                                                 GenomeId.index(crosslink[1].split(' ')[0]))
        elif crosslink[0].split(' ')[0] in CLASS2 and crosslink[1].split(' ')[0] in CLASS2:
            color = colors.linearlyInterpolatedColor(colors.purple, colors.red, 0, len(GenomeId),
                                                 GenomeId.index(crosslink[1].split(' ')[0]))
        else:
            color = colors.linearlyInterpolatedColor(colors.blue, colors.cyan, 0, len(GenomeId),
                                                 GenomeId.index(crosslink[1].split(' ')[0]))
        # color = list(colors.getAllNamedColors().keys())[GenomeId.index(crosslink[1].split(' ')[0]) * 17 + 17 % 163]
        F_x = set_X.add_feature(
            SeqFeature(FeatureLocation(int(crosslink[2]), int(crosslink[3]), strand=0)),
            color=color,
            border=color,
        )
        F_y = set_Y.add_feature(
            SeqFeature(FeatureLocation(int(crosslink[4]), int(crosslink[5]), strand=0)),
            color=color,
            border=color,
        )
        link_xy = CrossLink(F_x, F_y, color, color)
        gd_diagram.cross_track_links.append(link_xy)

    for record in records:
        gd_feature_set = feature_sets[record.id]

        # 矫正ori
        for feature in record.features:
            if feature.type == 'rep_origin':
                print(record.description + ' 的起始位点在:' + str(feature.location.start))
                record = record[feature.location.start:] + record[:feature.location.start]
                if record.features[0].strand == -1:
                    print('daole')
                    record = record.reverse_complement(id=True, name=True, description=True, features=True,
                                                       annotations=True, letter_annotations=True)
                break
        # 务必绘制反向互补序列时手动开启
        # record = record.reverse_complement(id=True, name=True, description=True, features=True,
        #                                    annotations=True, letter_annotations=True)

        print(record.description + ' 的起始位点已校正')
        # 画features
        i = 0
        if ref[record.id] != ['']:
            for feature in record.features:
                if feature.type != "gene":
                    continue
                color = list(colors.getAllNamedColors().keys())[len(feature) % 163]
                gd_feature_set.add_feature(feature, color=color, label=True, label_size=10, label_angle=90,
                                           sigil="ARROW", arrowshaft_height=1.0, arrowhead_length=0.1)
                i += 1
        elif ref[record.id] == ['']:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                color = list(colors.getAllNamedColors().keys())[len(feature) % 163]
                gd_feature_set.add_feature(feature, color=color, label=True, label_size=10, label_angle=90,
                                           sigil="ARROW", arrowshaft_height=1.0, arrowhead_length=0.2)
                i += 1
                # 用来手动添加重组位点
                # for pos in recombinations:
                #     if pos in record.features:
                #         gd_feature_set.add_feature(feature, color=color, label=True, label_size=10, label_angle=90,
                #                                    sigil="ARROW", arrowshaft_height=1.0, arrowhead_length=0.1)

    if not cross:
        # 用来绘制单一序列
        gd_diagram.draw(format="linear", pagesize='A4',
                        fragments=5, start=0, end=max_len, fragment_size=1)
        gd_diagram.write("T7.pdf", "PDF")
    else:
        # 用来绘制比对序列
        gd_diagram.draw(format="linear", pagesize=(10 * len(GenomeId) * cm, 120 * cm),
                        fragments=1, start=0, end=max_len, fragment_size=1)
        gd_diagram.write(output, "PDF")

    print("已输出为PDF")


data_write_csv('E:\\python_pycharm\data1.csv', new_data)
GenomeMap('E:\\python_pycharm\data1.csv', GenomeId)
