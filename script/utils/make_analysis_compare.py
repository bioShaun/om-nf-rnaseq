import pandas as pd
import fire
import itertools
from pathlib import Path


def rmats_sample_files(bam_dir, sample_group, out_dir, contrast=None):
    out_dir = Path(out_dir)
    bam_dir = Path(bam_dir).resolve()
    sample_df = pd.read_table(sample_group, header=None,
                              names=['sample_id'], index_col=1)
    contrast_obj = []
    if contrast is not None:
        contrast_df = pd.read_table(contrast, header=None,
                                    names=['g1', 'g2'])
        for each_idx in contrast_df.index:
            contrast_obj.append(list(contrast_df.loc[each_idx]))
    else:
        contrast_obj = itertools.combinations(
            sample_df.index.unique(), 2)

    for each_comp in contrast_obj:
        each_comp_name = f'{each_comp[0]}_vs_{each_comp[1]}'
        each_comp_path = out_dir / each_comp_name
        each_comp_path.mkdir(exist_ok=True, parents=True)
        for n, comp in enumerate(each_comp):
            comp_samples = sample_df.loc[comp, 'sample_id']
            bam_files = [str(bam_dir / f'{sample}.sorted.bam')
                         for sample in comp_samples]
            comp_file = each_comp_path / f'b{n+1}.txt'
            bam_file_string = ','.join(bam_files)
            with open(comp_file, 'w') as comp_bam_inf:
                comp_bam_inf.write(bam_file_string)


def generate_contrast(sample_group, contrast_file):
    sample_df = pd.read_table(sample_group, header=None,
                              names=['sample_id'], index_col=1)
    contrast_obj = itertools.combinations(
        sample_df.index.unique(), 2)
    with open(contrast_file, 'w') as contrast_inf:
        for each_com in contrast_obj:
            contrast_inf.write(f'{each_com[0]}\t{each_com[1]}\n')


if __name__ == '__main__':
    fire.Fire({'rmats_sample_files': rmats_sample_files,
               'generate_contrast': generate_contrast})
