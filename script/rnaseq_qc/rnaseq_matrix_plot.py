import fire
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO
from pathlib import PurePath

# set seaborn feature
sns.set()
sns.set_style("white")


PIE_COLOR1 = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
PIE_COLOR2 = ['#99cc66', '#ffff66', '#ff6666']
UTR_REGIONS = ['PCT_CODING_BASES', 'PCT_UTR_BASES',
               'PCT_INTRONIC_BASES', 'PCT_INTERGENIC_BASES']
NON_UTR_REGIONS = ['PCT_MRNA_BASES',
                   'PCT_INTRONIC_BASES', 'PCT_INTERGENIC_BASES']
STAR_COLS = ['Number of input reads',
             'Uniquely mapped reads number',
             'Number of reads mapped to multiple loci',
             'Number of reads mapped to too many loci']

COLOR_PAL = {
    'genome_region': PIE_COLOR1,
    'mapping': PIE_COLOR2,
}


def header2label(header):
    label = header.split('_')[1]
    if label == 'MRNA':
        label = 'EXON'
    return label


def pie_plot(labels, values, outfile_prefix, title,
             plot_content='genome_region'):
    val_num = len(values)
    fig1, ax1 = plt.subplots()
    colors = COLOR_PAL.get(plot_content)
    patches, texts, autotexts = ax1.pie(
        values, colors=colors[:val_num], labels=labels,
        autopct='%1.1f%%', startangle=90, pctdistance=0.8)

    for text in texts:
        text.set_color('grey')

    for autotext in autotexts:
        autotext.set_color('grey')

    ax1.axis('equal')
    ax1.set_title(title)
    plt.tight_layout()
    plt.savefig(f'{outfile_prefix}.pdf')
    plt.savefig(f'{outfile_prefix}.png')
    plt.close()


def genomic_origin_stats(matrix_file):
    genomic_origin_lines = ''
    flag = 0
    with open(matrix_file) as matrix_file_inf:
        for eachline in matrix_file_inf:
            if eachline.startswith('PF_BASES'):
                genomic_origin_lines += eachline
                flag = 1
            elif flag:
                genomic_origin_lines += eachline
                flag = 0
    genomic_origin_obj = StringIO(genomic_origin_lines)
    genomic_origin_df = pd.read_csv(genomic_origin_obj, sep='\t')
    if genomic_origin_df.loc[0].PCT_UTR_BASES > 0:
        labels = UTR_REGIONS
    else:
        labels = NON_UTR_REGIONS
    values = genomic_origin_df.loc[0, labels].values
    labels = list(map(header2label, labels))
    return labels, values


def rnaseq_matrix_plot(matrix_file):
    matrix_file = PurePath(matrix_file)
    matrix_prefix = matrix_file.with_suffix('')
    sample_name = matrix_prefix.stem
    # plot reads distribution in genome region
    pie_plot_prefix = matrix_prefix.with_suffix('.genome_region')
    labels, values = genomic_origin_stats(matrix_file)
    pie_plot(labels, values, pie_plot_prefix, sample_name)
    # plot reads cov
    cov_df = pd.read_csv(matrix_file,
                         sep='\t', skiprows=list(range(10)))
    cov_df.columns = ['position', 'coverage']
    ax = sns.lineplot(x="position", y="coverage",
                      markers=True, dashes=False, data=cov_df)
    ax.set_title(sample_name)
    fig = ax.get_figure()
    outprefix = matrix_prefix.with_suffix('.gene_coverage')
    fig.savefig(f'{outprefix}.png')
    fig.savefig(f'{outprefix}.pdf')


def mapping_plot(star_log, sample_name):
    star_df = pd.read_csv(star_log, header=None, sep='|', index_col=0)
    star_df.index = [each.strip() for each in star_df.index]
    star_df = star_df.loc[STAR_COLS]
    star_df.loc[:, 1] = [each.strip() for each in star_df.loc[:, 1]]
    star_stats = list(map(int, star_df.loc[STAR_COLS, 1].values))
    areads, unique_mapped, multi_mapped1, multi_mapped2 = star_stats
    multi_mapped = multi_mapped1 + multi_mapped2
    unmapped = areads - unique_mapped - multi_mapped
    labels = ['Unique mapped', 'Multiple mapped', 'Unmapped']
    values = [unique_mapped, multi_mapped, unmapped]
    outfile = PurePath(star_log).parent / f'{sample_name}.mapping_rate'
    pie_plot(labels, values, outfile, sample_name, plot_content='mapping')


if __name__ == '__main__':
    fire.Fire({'rnaseq_matrix': rnaseq_matrix_plot,
               'mapping': mapping_plot})
