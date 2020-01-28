import fire
import gtfparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

ANN_FEATURE_NAMES = ('protein_coding', 'annotated_lncRNA')
PLOT_ORDER = ('novel_lncRNA', 'annotated_lncRNA', 'protein_coding')
PIE_COLOR_MAP = {
    'novel_lncRNA': '#e74c3c',
    'TUCP': '#3498db',
    'annotated_lncRNA': '#34495e',
    'protein_coding': '#2ecc71',
    'other_ncRNA': '#95a5a6',
    'pseudogene': '#9b59b6',
}


# set seaborn feature
sns.set()
sns.set_style("white")


def pie_plot(labels, values, outfile_prefix, title, color_pal):
    val_num = len(values)
    fig1, ax1 = plt.subplots()
    patches, texts = ax1.pie(
        values, colors=color_pal[:val_num],
        startangle=90, pctdistance=0.8)

    plt.legend(patches, labels, loc="best")
    ax1.axis('equal')
    ax1.set_title(title)
    plt.tight_layout()
    plt.savefig(f'{outfile_prefix}.pdf')
    plt.savefig(f'{outfile_prefix}.png')
    plt.close()


def cal_plot_limit(df, label, plot_cut=0.8):
    plot_limit = df.loc[:, label].quantile(plot_cut)
    plot_limit = int(np.round(plot_limit, int(-np.log10(plot_limit))))
    return plot_limit


def density_plot(df, label, outdir, plot_cut=0.8):
    plot_limit = cal_plot_limit(df, label, plot_cut)
    for gt in PLOT_ORDER:
        # Subset to the airline
        subset = df[df.gene_biotype == gt]
        subset = subset[subset[label] <= plot_limit]
        # Draw the density plot
        ax = sns.distplot(subset[label], hist=False, kde=True,
                          kde_kws={'linewidth': 1.5, 'shade': True},
                          label=gt)
    xlab = label.replace('_', ' ')
    ax.set_xlabel(xlab)
    fig = ax.get_figure()
    fig.set_size_inches(8, 6)
    out_prefix = outdir / f'lncRNA_feature_{label.lower()}'
    fig.savefig(f'{out_prefix}.png')
    fig.savefig(f'{out_prefix}.pdf')
    fig.clf()


def barplot(df, label, outdir, plot_cut=0.8):
    plot_df = df.copy()
    plot_limit = cal_plot_limit(df, label, plot_cut)
    mask_idx = plot_df[plot_df[label] >= plot_limit].index
    plot_df.loc[mask_idx, label] = f'≥{plot_limit}'
    plot_stats = plot_df.groupby([label, 'gene_biotype']).size()
    plot_portion = plot_stats / plot_stats.sum(level=1)
    plot_portion.name = 'Proportion'
    plot_portion = pd.DataFrame(plot_portion).reset_index()
    ax = sns.barplot(x=label, y="Proportion",
                     hue="gene_biotype",
                     edgecolor=".6", data=plot_portion,
                     hue_order=PLOT_ORDER)
    ax.legend(loc='upper right')
    xlab = label.replace('_', ' ')
    ax.set_xlabel(xlab)
    fig = ax.get_figure()
    fig.set_size_inches(8, 6)
    out_prefix = outdir / f'lncRNA_feature_{label.lower()}'
    fig.savefig(f'{out_prefix}.png')
    fig.savefig(f'{out_prefix}.pdf')
    fig.clf()


def test_plot(plot_file, outdir):
    outdir = Path(outdir)
    gene_inf_df = pd.read_csv(plot_file)
    # barplot_stats(gene_inf_df, 'Isoform_number', outdir, plot_cut=0.99)
    # barplot_stats(gene_inf_df, 'Exon_number', outdir)
    # density_plot(gene_inf_df, 'Length', outdir)
    # gene number plot
    novel_gene_df = gene_inf_df[gene_inf_df.gene_biotype.isin(
        ('TUCP', 'novel_lncRNA', 'lncRNA'))]
    novel_gene_pre = outdir / 'Novel_gene_number'
    gene_number_plot(novel_gene_df, novel_gene_pre)
    all_gene_pre = outdir / 'All_gene_number'
    gene_number_plot(gene_inf_df, all_gene_pre)


def gene_number_plot(df, out_prefix):
    type_stats = df.gene_biotype.value_counts()
    type_stats_por = type_stats / len(df)
    labels = []
    colors = []
    for n, idx_i in enumerate(type_stats.index):
        n_count = type_stats.values[n]
        n_por = type_stats_por.values[n]
        labels.append(f'{idx_i} ({n_count:,}, {n_por:.1%})')
        colors.append(PIE_COLOR_MAP.get(idx_i))
    pie_plot(labels, type_stats.values, out_prefix,
             title="", color_pal=colors)


def gene_locus(gtf_df):
    chrom = gtf_df.groupby(['gene_id', 'gene_biotype'])['seqname'].first()
    start = gtf_df.groupby(['gene_id', 'gene_biotype'])['start'].min()
    end = gtf_df.groupby(['gene_id', 'gene_biotype'])['end'].max()
    strand = gtf_df.groupby(['gene_id', 'gene_biotype'])['strand'].first()
    locus = pd.concat([chrom, start, end, strand], axis=1)
    locus = locus.astype('str')
    locus.loc[:, 'Locus'] = locus.seqname.str.cat(locus.start, sep=':')
    locus.loc[:, 'Locus'] = locus.Locus.str.cat(locus.end, sep='-')
    locus.loc[:, 'Locus'] = locus.Locus.str.cat(locus.strand, sep='|')
    return locus.Locus


def lnc_feature(novel_lnc, gtf_split_dir, outdir, gene_ann=None):
    gtf_df_list = []
    gtf_df_list.append(gtfparse.read_gtf(novel_lnc))
    split_gtfs = Path(gtf_split_dir).glob('*gtf')
    for gtf in split_gtfs:
        gtf_df_list.append(gtfparse.read_gtf(gtf))
    merge_gtf_df = pd.concat(gtf_df_list)
    exon_df = merge_gtf_df[merge_gtf_df.feature == 'exon']
    # summarize tr feature
    tr_exon_num = exon_df.groupby(['transcript_id']).size()
    tr_exon_num.name = 'Exon_number'
    exon_df.loc[:, 'Length'] = exon_df.end - exon_df.start + 1
    tr_len_df = exon_df.groupby(['transcript_id'])['Length'].sum()
    tr_inf_df = pd.concat([tr_len_df, tr_exon_num], axis=1)
    tr_inf_df = tr_inf_df.reset_index()
    # summarize gene feature
    gene_locus_df = gene_locus(exon_df)
    gene_tr_df = exon_df.loc[:, [
        'gene_id', 'transcript_id',
        'gene_biotype']].drop_duplicates()
    gene_tr_df = gene_tr_df.merge(tr_inf_df, sort=False)
    gene_isoform = gene_tr_df.groupby(['gene_id', 'gene_biotype']).size()
    gene_isoform.name = 'Isoform_number'
    gene_length = gene_tr_df.groupby(['gene_id', 'gene_biotype'])[
        'Length'].max()
    gene_exon = gene_tr_df.groupby(['gene_id', 'gene_biotype'])[
        'Exon_number'].max()
    gene_inf = pd.concat([gene_locus_df, gene_isoform,
                          gene_length, gene_exon], axis=1)
    gene_inf_df = gene_inf.reset_index()
    if gene_ann is not None:
        gene_ann_df = pd.read_csv(gene_ann, sep='\t')
        gene_inf_df = gene_inf_df.merge(gene_ann_df, how='left')
    outdir = Path(outdir)
    gene_feature_file = outdir / 'Gene_feature.csv'
    gene_inf_df.to_csv(gene_feature_file, na_rep='--', index=False)
    # feature plot
    feature_plot_df = gene_inf_df[gene_inf_df.gene_biotype.isin(PLOT_ORDER)]
    barplot(feature_plot_df, 'Isoform_number', outdir, plot_cut=0.99)
    barplot(feature_plot_df, 'Exon_number', outdir)
    density_plot(feature_plot_df, 'Length', outdir)
    # gene number plot
    novel_gene_df = gene_inf_df[gene_inf_df.gene_biotype.isin(
        ('TUCP', 'novel_lncRNA', 'lncRNA'))]
    novel_gene_pre = outdir / 'Novel_gene_number'
    gene_number_plot(novel_gene_df, novel_gene_pre)
    all_gene_pre = outdir / 'All_gene_number'
    gene_number_plot(gene_inf_df, all_gene_pre)


if __name__ == '__main__':
    fire.Fire(lnc_feature)
