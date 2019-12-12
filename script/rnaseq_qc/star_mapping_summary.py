import fire
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

sns.set()
sns.set_style("white")


def read_star_mapping_log(star_log_file):
    star_df = pd.read_table(star_log_file, header=None, sep='|', index_col=0)
    star_df = star_df.dropna()
    star_df = star_df.iloc[4:]
    star_df.loc[:, 1] = [each.strip() for each in star_df.loc[:, 1]]
    return star_df


def mapping_plot(star_df, out_prefix):
    plot_df = star_df.reset_index()

    def percent2num(per_num):
        return float(per_num.rstrip('%'))

    plot_df.loc[:, 'Unique_mapped'] = plot_df.loc[
        :, 'Uniquely mapped reads %'].map(percent2num)
    plot_df.loc[:, 'multi_mapped1'] = plot_df.loc[
        :, '% of reads mapped to multiple loci'].map(percent2num)
    plot_df.loc[:, 'multi_mapped2'] = plot_df.loc[
        :, '% of reads mapped to too many loci'].map(percent2num)
    plot_df.loc[:, 'Multiple_mapped'] = plot_df.multi_mapped1 + \
        plot_df.multi_mapped2
    m_plot_df = plot_df.melt(id_vars=['Sample'],
                             value_vars=['Unique_mapped',
                                         'Multiple_mapped'],
                             value_name='percent_of_reads',
                             var_name='mapping_stats')
    ax = sns.lineplot(x="Sample", dashes=False,
                      y="percent_of_reads", hue="mapping_stats",
                      style="mapping_stats", markers=["o", "o"],
                      data=m_plot_df)
    fig = ax.get_figure()
    fig.savefig(f'{out_prefix}.png')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,
                       horizontalalignment='right')
    ax.set_xlabel('')
    ax.set_ylabel('Percent of Reads')
    sample_num = len(ax.get_xticklabels())
    plot_height = 8 + sample_num / 10
    plot_width = 6 + sample_num / 5
    fig = ax.get_figure()
    fig.set_size_inches(plot_width, plot_height)
    fig.subplots_adjust(bottom=0.2)
    fig.savefig(f'{out_prefix}.png')
    fig.savefig(f'{out_prefix}.pdf')


def star_mapping_summary(star_log_dir, out_dir):
    star_log_file_list = list(Path(star_log_dir).glob('*.Log.final.out'))
    star_log_df_list = list(map(read_star_mapping_log, star_log_file_list))
    star_log_df = pd.concat(star_log_df_list, axis=1)
    star_log_df.columns = [each.name.replace(
        '.Log.final.out', '') for each in star_log_file_list]
    star_log_out_df = star_log_df.T
    star_log_out_df.index.name = 'Sample'
    star_log_out_df.columns = [each.strip()
                               for each in star_log_out_df.columns]

    # output
    star_summary_file = Path(out_dir) / 'mapping_summary.csv'
    star_log_out_df.to_csv(star_summary_file)
    mapping_plot(star_log_out_df, star_summary_file.stem)


if __name__ == '__main__':
    fire.Fire(star_mapping_summary)
