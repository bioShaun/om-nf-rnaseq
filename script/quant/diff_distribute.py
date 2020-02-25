import fire
import delegator
import pandas as pd
from pathlib import PurePath
from pybedtools import BedTool


R_BIN = '/public/software/R/R-3.5.1/executable/bin/'
PLOT_R = PurePath(__file__).parent / 'diff_distribution.R'


def diff_dis(diff_file, chr_size, chr_window, pval=0.05, logfc=1, lnc=False):
    diff_file = PurePath(diff_file)
    compare_name = diff_file.name.replace('.edgeR.DE_results.txt', '')
    outdir = diff_file.parent
    diff_df = pd.read_csv(diff_file, sep='\t')
    diff_df = diff_df[diff_df.Locus != "--"]
    diff_df.loc[:, "chrom"] = [each.split(':')[0] for each in diff_df.Locus]
    diff_df.loc[:, 'end'] = [each.split(':')[1].split(
        '|')[0].split('-')[1] for each in diff_df.Locus]
    diff_df.loc[:, 'start'] = [int(each.split(':')[1].split(
        '|')[0].split('-')[0]) - 1 for each in diff_df.Locus]
    mask1 = diff_df.FDR <= pval
    mask2 = diff_df.logFC.abs() >= logfc
    f_diff_df = diff_df[mask1 & mask2]
    plot_by_type = ''
    if lnc:
        mask = f_diff_df.gene_biotype.str.contains('lncRNA')
        f_diff_df.loc[f_diff_df[mask].index, 'gene_biotype'] = 'lncRNA'
        mask1 = f_diff_df.gene_biotype == 'protein_coding'
        mask2 = f_diff_df.gene_biotype == 'lncRNA'
        plot_by_type = '--lnc'
    else:
        f_diff_df.loc[:, 'gene_biotype'] = 'protein_coding'
    pcg_lnc_diff_df = f_diff_df[mask1 | mask2]
    if pcg_lnc_diff_df.empty:
        return 'Zero Diff gene Found'
    pcg_lnc_diff_df.loc[:, 'regulation'] = 'UP'
    mask = pcg_lnc_diff_df.logFC < 0
    if sum(mask):
        pcg_lnc_diff_df.loc[pcg_lnc_diff_df[mask].index, 'regulation'] = 'DOWN'
    gene_bed_file = outdir / 'diff.gene.bed'
    pcg_lnc_diff_df.to_csv(
        gene_bed_file, sep='\t',
        header=False, columns=[
            'chrom', 'start', 'end',
            'gene_biotype', 'regulation'], index=False)
    gene_bed = BedTool(str(gene_bed_file))
    window_bed = BedTool(chr_window)
    intersect_obj = window_bed.intersect(gene_bed, wo=True)
    intersect_df = intersect_obj.to_dataframe()
    intersect_count = intersect_df.groupby(
        ['chrom', 'start', 'end', 'thickStart', 'thickEnd']).size()
    intersect_count.name = 'gene_count'
    intersect_count_df = pd.DataFrame(intersect_count)
    intersect_count_file = outdir / 'diff_gene_count.csv'
    intersect_count_df.to_csv(intersect_count_file)
    plot_cmd = (f'{R_BIN}/Rscript {PLOT_R} '
                f'--diff_file {intersect_count_file} '
                f'--chrom_size {chr_size} '
                f'--out_dir {outdir} '
                f'--compare {compare_name} {plot_by_type}')
    # clean_cmd = f'rm {intersect_count_file} {gene_bed_file}'
    delegator.run(plot_cmd)


if __name__ == '__main__':
    fire.Fire(diff_dis)
