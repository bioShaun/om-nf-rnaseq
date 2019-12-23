import fire
import numpy as np
import pandas as pd
from pathlib import Path


def de_gene_num(de_gene_file):
    if de_gene_file.is_file():
        de_gene_df = pd.read_csv(de_gene_file, header=None)
        return len(de_gene_df)
    else:
        return 0


def de_number(diff_dir, group_inf, contrasts, outdir):
    outdir = Path(outdir)
    group_df = pd.read_csv(group_inf,
                           index_col=0,
                           sep='\t',
                           header=None,
                           names=['group_id'])
    groups = group_df.group_id.unique()
    diff_num_df = pd.DataFrame([[np.nan for each in groups]
                                for each in groups],
                               index=groups,
                               columns=groups)
    compare_num_df = pd.DataFrame([], columns=['compare', 'up', 'down'])
    diff_dir = Path(diff_dir)
    for n, diff_i in enumerate(diff_dir.iterdir()):
        compare_i = diff_i.name
        upname, dname = compare_i.split('_vs_')
        up_gene = diff_i / \
            f'{compare_i}.{upname}-UP.edgeR.DE_results.diffgenes.txt'
        down_gene = diff_i / \
            f'{compare_i}.{dname}-UP.edgeR.DE_results.diffgenes.txt'
        up_gene_num = de_gene_num(up_gene)
        down_gene_num = de_gene_num(down_gene)
        diff_num_df.loc[upname, dname] = up_gene_num
        diff_num_df.loc[dname, upname] = down_gene_num
        compare_num_df.loc[n] = [compare_i, up_gene_num, down_gene_num]
    if contrasts:
        for group_i in groups:
            diff_num_df.loc[group_i, group_i] = 0
        diff_num_df = diff_num_df.astype('int')
        outfile = outdir / 'diff.genes.matrix.csv'
        diff_num_df.to_csv(outfile)
    else:
        outfile = outdir / 'diff.genes.compare.csv'
        compare_num_df.to_csv(outfile, index=False)


if __name__ == '__main__':
    fire.Fire(de_number)
