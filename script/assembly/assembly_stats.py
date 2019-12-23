import fire
import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

# set seaborn feature
sns.set()
sns.set_style("white")

ANNOTATED_CLS = ['=', 'c', 'k', 'm', 'n', 'j', 'e', 'o']
ERR_CLS = 'y'
INTERGENIC_CODE = ['u', 'e']

COLOR_PAL = ['#c3f584', '#ffd271']


def pie_plot(labels, values, outfile_prefix, title):
    val_num = len(values)
    fig1, ax1 = plt.subplots()
    colors = COLOR_PAL
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


def line_plot(plt_df, x, y, cat, out_prefix, xlab='', ylab='', ylim=None):
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
    if ylim:
        ax.set_ylim(ylim)
    sample_num = len(ax.get_xticklabels())
    plot_height = 8 + sample_num / 10
    plot_width = 6 + sample_num / 5
    fig = ax.get_figure()
    fig.set_size_inches(plot_width, plot_height)
    fig.subplots_adjust(bottom=0.2)
    fig.savefig(f'{out_prefix}.png')
    fig.savefig(f'{out_prefix}.pdf')
    fig.clf()


def riboz_cls_code_des(cls_code):
    if cls_code in ANNOTATED_CLS:
        return 'Annotated'
    else:
        return 'Unannotated'


def polya_cls_code_des(cls_code):
    if cls_code in INTERGENIC_CODE:
        return 'Unannotated'
    else:
        return 'Annotated'


def gffcompare_stats(tmap, outdir, stranded=False, plot=True):
    tmap_df = pd.read_csv(tmap, sep='\t')
    tmap_df = tmap_df[tmap_df.class_code != ERR_CLS]
    if stranded:
        tmap_df.loc[:, 'class_code'] = tmap_df.class_code.map(
            riboz_cls_code_des)
    else:
        tmap_df.loc[:, 'class_code'] = tmap_df.class_code.map(
            polya_cls_code_des)
    ann_stats = tmap_df.class_code.value_counts()
    if plot:
        plot_pre = Path(outdir) / 'assembly_summary'
        labels = []
        for n, idx_i in enumerate(ann_stats.index):
            labels.append(f'{idx_i} ({ann_stats.values[n]:,})')
        pie_plot(labels, ann_stats.values, plot_pre, title="")
    else:
        filename = tmap.stem.replace('cmp2ref.', '').replace('.gtf', '')
        ann_stats.name = filename
        return ann_stats


def gffcompare_stats_grp(tmap_dir, out_dir, stranded=False):
    tmap_files = Path(tmap_dir).glob('*.tmap')
    tmap_df_list = []
    for tmap_i in tmap_files:
        tmap_i_df = gffcompare_stats(
            tmap_i, None, stranded=stranded, plot=False)
        tmap_df_list.append(tmap_i_df)
    tmap_df = pd.concat(tmap_df_list, axis=1).T
    tmap_df.index.name = 'Sample_id'
    tmap_df = tmap_df.reset_index()
    plot_df = tmap_df.melt(id_vars='Sample_id',
                           value_name='Transcripts', var_name='Type')
    out_pre = Path(out_dir) / 'assembly_summary_by_sample'
    line_plot(plot_df, 'Sample_id', 'Transcripts', 'Type', out_pre,
              ylab='Number of transcripts')


if __name__ == '__main__':
    fire.Fire({'stats_combined': gffcompare_stats,
               'stats_all': gffcompare_stats_grp})
