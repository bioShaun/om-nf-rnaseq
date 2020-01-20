
import fire
import pandas as pd
from gtfparse import read_gtf
from pathlib import PurePath, Path

LABLE2TYPE_DICT = {
    'noncoding': 'lncRNA',
    'coding': 'TUCP',
}

LABLE2TYPE_DICT2 = {
    'noncoding': 'novel_lncRNA',
    'coding': 'TUCP',
}


def label2genetype(unique_label):
    if len(unique_label) == 2:
        return 'TUCP'
    else:
        return unique_label[0]


def dfline2gtfline(dfline):
    basic_inf = dfline[:8]
    basic_inf.fillna('.', inplace=True)
    basic_inf.frame = '.'
    basic_inf_list = [str(each) for each in basic_inf]
    basic_inf_line = '\t'.join(basic_inf_list)
    attr_inf = dfline[8:]
    attr_inf_list = []
    for key, val in attr_inf.items():
        if val:
            attr_inf_list.append(f'{key} "{val}";')
    attr_inf_line = ' '.join(attr_inf_list)
    return f'{basic_inf_line}\t{attr_inf_line}\n'


def cpc2lnc(cpc, gtf, split_gtf_dir, outfile):
    gtf_df = read_gtf(gtf)
    tr_df = gtf_df[gtf_df.feature == 'transcript']
    tr2gene_df = tr_df.loc[:, ['transcript_id', 'gene_id']]
    cpc_df = pd.read_csv(cpc, sep='\t')
    gene_cpc_df = cpc_df.merge(
        tr2gene_df, left_on='#ID', right_on='transcript_id')
    ann_gtf = Path(split_gtf_dir) / 'annotated_lncRNA.gtf'
    if ann_gtf.is_file():
        gene_cpc_df.loc[:, 'transcript_biotype'] = gene_cpc_df.label.map(
            LABLE2TYPE_DICT2)
    else:
        gene_cpc_df.loc[:, 'transcript_biotype'] = gene_cpc_df.label.map(
            LABLE2TYPE_DICT)
    gene_biotype = gene_cpc_df.groupby(
        ['gene_id'])['transcript_biotype'].unique()
    gene_biotype = gene_biotype.map(label2genetype)
    gene_biotype.name = 'gene_biotype'
    gene_biotype = pd.DataFrame(gene_biotype)
    gene_cpc_df = gene_cpc_df.merge(
        gene_biotype, left_on='gene_id', right_on='gene_id')
    type_df = gene_cpc_df.loc[:, ['transcript_id',
                                  'gene_id',
                                  'transcript_biotype',
                                  'gene_biotype']]
    type_gtf_df = gtf_df.merge(type_df)
    with open(outfile, 'w') as type_gtf_inf:
        for idx in type_gtf_df.index:
            outline = dfline2gtfline(type_gtf_df.loc[idx])
            type_gtf_inf.write(outline)


if __name__ == '__main__':
    fire.Fire(cpc2lnc)
