import click
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO
from pathlib import Path


# set seaborn feature
sns.set()
sns.set_style("white")

UTR_REGIONS = ['PCT_CODING_BASES', 'PCT_UTR_BASES',
               'PCT_INTRONIC_BASES', 'PCT_INTERGENIC_BASES']
NON_UTR_REGIONS = ['PCT_MRNA_BASES',
                   'PCT_INTRONIC_BASES', 'PCT_INTERGENIC_BASES']


RNASEQ_METIRCS_PATTERN = '.RnaSeqMetrics.txt'
RNASEQ_METIRCS_HEADER = 'PF_BASES'

IS_METIRCS_PATTERN = '.insert_size_metrics.txt'
IS_METIRCS_HEADER = 'MEDIAN_INSERT_SIZE'


def get_metrix_inf(m_file, file_suffix, header_label):
    sample_name = m_file.name.replace(
        file_suffix, '')
    matrix_inf = ''
    with open(m_file) as file_inf:
        flag = 0
        for eachline in file_inf:
            if header_label in eachline:
                flag = 1
                matrix_inf += eachline
            elif flag:
                matrix_inf += eachline
                break
    genomic_origin_obj = StringIO(matrix_inf)
    genomic_origin_df = pd.read_csv(genomic_origin_obj, sep='\t')
    genomic_origin_df.index = [sample_name]
    return genomic_origin_df


def header2label(header):
    label = header.split('_')[1]
    if label == 'MRNA':
        label = 'EXON'
    return label


def line_plot(plt_df, x, y, cat, out_prefix, xlab='', ylab=''):
    cat_num = len(plt_df.loc[:, cat].unique())
    markers = ["o"] * cat_num
    ax = sns.lineplot(x=x, dashes=False,
                      y=y, hue=cat,
                      style=cat, markers=markers,
                      data=plt_df)
    fig = ax.get_figure()
    fig.savefig(f'{out_prefix}.png')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,
                       horizontalalignment='right')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    sample_num = len(ax.get_xticklabels())
    plot_height = 8 + sample_num / 10
    plot_width = 6 + sample_num / 5
    fig = ax.get_figure()
    fig.set_size_inches(plot_width, plot_height)
    fig.subplots_adjust(bottom=0.2)
    fig.savefig(f'{out_prefix}.png')
    fig.savefig(f'{out_prefix}.pdf')
    fig.clf()


def rnaseq_metrics(analysis_dir, out_dir):
    summary_file = Path(out_dir) / 'RnaSeqMetrics_summary.csv'
    metrics_files = Path(analysis_dir).glob(f'*{RNASEQ_METIRCS_PATTERN}')
    metrics_df_list = [get_metrix_inf(
        each, RNASEQ_METIRCS_PATTERN, RNASEQ_METIRCS_HEADER)
        for each in metrics_files]
    metrics_df = pd.concat(metrics_df_list)
    metrics_df.index.name = 'Sample'
    metrics_df.drop(['SAMPLE', 'LIBRARY', 'READ_GROUP'],
                    inplace=True,
                    axis=1)
    metrics_df.to_csv(summary_file)
    metrics_df = metrics_df.reset_index()
    # genome region plot
    if metrics_df.PCT_UTR_BASES.min() > 0:
        region_list = UTR_REGIONS
    else:
        region_list = NON_UTR_REGIONS
    genome_region_df = metrics_df.melt(id_vars=['Sample'],
                                       value_vars=region_list,
                                       value_name='percent_of_reads',
                                       var_name='genome_region')
    genome_region_df.loc[:, 'genome_region'] = genome_region_df.genome_region.map(
        header2label)
    genome_region_pre = Path(out_dir) / 'overall_genome_region'
    line_plot(genome_region_df, 'Sample', 'percent_of_reads',
              'genome_region', genome_region_pre,
              ylab='Percent of Reads')

    # gene coverage plot
    cov_df = metrics_df.melt(id_vars=['Sample'],
                             value_vars=['MEDIAN_5PRIME_TO_3PRIME_BIAS'],
                             value_name='value',
                             var_name='Statistics')
    cov_pre = Path(out_dir) / 'overall_gene_coverage'
    line_plot(cov_df, 'Sample', 'value',
              'Statistics', cov_pre,
              ylab="Ratio of coverage at 5' to 3'")


def is_metrics(analysis_dir, out_dir):
    summary_file = Path(out_dir) / 'InsertSizeMetrics_summary.csv'
    metrics_files = Path(analysis_dir).glob(f'*{IS_METIRCS_PATTERN}')
    metrics_df_list = [get_metrix_inf(
        each, IS_METIRCS_PATTERN, IS_METIRCS_HEADER)
        for each in metrics_files]
    metrics_df = pd.concat(metrics_df_list)
    metrics_df.index.name = 'Sample'
    metrics_df.to_csv(summary_file)
    metrics_df = metrics_df.reset_index()
    is_df = metrics_df.melt(id_vars=['Sample'],
                            value_vars=['MEDIAN_INSERT_SIZE'],
                            value_name='insert_size',
                            var_name='Statistics')
    is_pre = Path(out_dir) / 'overall_insert_size'
    line_plot(is_df, 'Sample', 'insert_size',
              'Statistics', is_pre,
              ylab="Insert Size")


@click.command()
@click.option('-d', '--analysis_dir', type=click.Path(exists=True),
              help='rnaseq matrix directory.', required=True)
@click.option('-o', '--out_dir', type=click.Path(exists=True),
              help='rnaseq matrix summary directory.', required=True)
@click.option('-m', '--metrics', type=click.Choice(['rnaseq', 'is']),
              help='metrics type to extract info', required=True)
def main(analysis_dir, out_dir, metrics):
    if metrics == 'rnaseq':
        rnaseq_metrics(analysis_dir, out_dir)
    elif metrics == 'is':
        is_metrics(analysis_dir, out_dir)
    else:
        click.secho(f'Unknown metrics {metrics}', fg='red')


if __name__ == '__main__':
    main()
