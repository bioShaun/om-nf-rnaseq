import fire
import json
import glob
from pandas import DataFrame
import pandas as pd
from pathlib import Path, PurePath


OUT_COL = [
    'raw_total_bases',
    'raw_total_reads',
    'total_bases',
    'total_reads',
    'clean_bases_rate',
    'clean_reads_rate',
    'q30_rate',
    'gc_content'
]

OUT_COL_MAP = {
    'raw_total_reads': 'Raw_reads',
    'raw_total_bases': 'Raw_bases',
    'total_reads': 'Clean_reads',
    'total_bases': 'Clean_bases',
    'clean_reads_rate': 'Clean_reads_rate',
    'clean_bases_rate': 'Clean_bases_rate',
    'q30_rate': 'Q30',
    'gc_content': 'GC'
}


def json2gc(fastp_json, reads_num=1):
    gc_label = f'read{reads_num}_after_filtering'
    gc_df = pd.DataFrame(fastp_json[gc_label]['content_curves'])
    gc_df.index.name = '#Base'
    gc_df = gc_df * 100
    gc_df = gc_df.reset_index()
    gc_df.drop('GC', axis=1, inplace=True)
    if reads_num == 1:
        off_set = 1
    else:
        off_set = 151
    gc_df.loc[:, '#Base'] = gc_df.loc[:, '#Base'] + off_set
    return gc_df


def json2quality(fastp_json, reads_num=1):
    data_label = f'read{reads_num}_after_filtering'
    quality_df = pd.DataFrame(
        fastp_json[data_label]['quality_curves']['mean'], columns=['quality'])
    quality_df.index.name = '#Base'
    quality_df = quality_df.reset_index()
    if reads_num == 1:
        off_set = 1
    else:
        off_set = 151
    quality_df.loc[:, '#Base'] = quality_df.loc[:, '#Base'] + off_set
    return quality_df


def extract_fastp_json(fastp_dir, outdir):
    outdir = Path(outdir)
    gc_dir = outdir / 'reads_gc'
    qual_dir = outdir / 'reads_quality'
    filter_dir = outdir / 'reads_filter'
    filter_dir.mkdir(parents=True, exist_ok=True)
    qual_dir.mkdir(parents=True, exist_ok=True)
    gc_dir.mkdir(parents=True, exist_ok=True)
    json_files = glob.glob(f'{fastp_dir}/*.json')
    fastp_dfs = []
    raw_dfs = []
    for each_json_file in json_files:
        with open(each_json_file) as json_inf:
            json_obj = json.load(json_inf)
            file_name = PurePath(each_json_file).stem
            each_df = pd.DataFrame(
                list(json_obj['summary']['after_filtering'].values()),
                index=list(json_obj['summary']['after_filtering'].keys()),
                columns=[file_name])
            fastp_dfs.append(each_df)
            raw_df = pd.DataFrame(
                list(json_obj['summary']['before_filtering'].values()),
                index=list(json_obj['summary']['before_filtering'].keys()),
                columns=[file_name]
            )
            raw_dfs.append(raw_df)
            read1_gc_df = json2gc(json_obj, 1)
            read2_gc_df = json2gc(json_obj, 2)
            gc_df = pd.concat([read1_gc_df, read2_gc_df])
            gc_file = gc_dir / f'{file_name}.gc.txt'
            gc_df.to_csv(gc_file, sep='\t', index=False)
            read1_qual_df = json2quality(json_obj, 1)
            read2_qual_df = json2quality(json_obj, 2)
            qual_df = pd.concat([read1_qual_df, read2_qual_df])
            qual_file = qual_dir / f'{file_name}.reads_quality.txt'
            qual_df.to_csv(qual_file, sep='\t', index=False)
            filter_file = filter_dir / f'{file_name}.filter.txt'
            with open(filter_file, 'w') as filter_inf:
                for each_key, each_val in json_obj['filtering_result'].items():
                    filter_inf.write(f'{each_key}\t{each_val}\n')
    # main summary
    fastp_clean_df = pd.concat(fastp_dfs, axis=1).T
    fastp_raw_df = pd.concat(raw_dfs, axis=1).T
    fastp_raw_df.columns = ['raw_{}'.format(
        each) for each in fastp_raw_df.columns]
    merged_df = fastp_raw_df.merge(fastp_clean_df,
                                   left_index=True,
                                   right_index=True)
    merged_df.loc[:, 'clean_bases_rate'] = merged_df.total_bases / \
        merged_df.raw_total_bases
    merged_df.loc[:, 'clean_reads_rate'] = merged_df.total_reads / \
        merged_df.raw_total_reads
    out_merged_df = merged_df.loc[:, OUT_COL]
    out_merged_df.columns = out_merged_df.columns.map(OUT_COL_MAP)
    out_merged_df.index.name = 'Sample_id'
    out_merged_df.loc[:, 'Raw_reads'] = out_merged_df.Raw_reads.astype('int')
    out_merged_df.loc[:, 'Raw_bases'] = out_merged_df.Raw_bases.astype('int')
    out_merged_df.loc[
        :, 'Clean_reads'] = out_merged_df.Clean_reads.astype('int')
    out_merged_df.loc[
        :, 'Clean_bases'] = out_merged_df.Clean_bases.astype('int')
    fastp_clean_file = outdir / 'data.summary.csv'
    out_merged_df.to_csv(fastp_clean_file,
                         float_format='%.3f')


if __name__ == '__main__':
    fire.Fire(extract_fastp_json)
