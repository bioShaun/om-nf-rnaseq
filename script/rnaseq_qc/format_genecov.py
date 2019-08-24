import fire
import pandas as pd
from pathlib import PurePath


def genecov4multiqc(genecov):
    genecov = PurePath(genecov)
    genecov_df = pd.read_csv(genecov, sep='\t')
    genecov_df.loc[:, 'percentile'] = genecov_df.percentile + 1
    genecov_df = genecov_df.set_index('percentile')
    genecov_df.columns = [genecov.name.replace('.geneBodyCoverage.txt', '')]
    t_genecov_df = genecov_df.T
    t_genecov_df.index.name = 'Percentile'
    t_genecov_df.to_csv(genecov, sep='\t')


if __name__ == '__main__':
    fire.Fire(genecov4multiqc)
