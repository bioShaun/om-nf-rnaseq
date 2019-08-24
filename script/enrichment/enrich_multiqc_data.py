import re
import fire
import json
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict


ENRICH_COL = {
    'kegg': ['#Term', 'Corrected P-Value'],
    'go': ['term', 'qvalue'],
}


def top_enriched_items(enrich_dir, enrich_type, out_dir, number=20):
    file_suffix = f'.ALL.{enrich_type}.enrichment.txt'
    pattern = re.compile('(\S+)%s' % file_suffix)
    enrich_dir = Path(enrich_dir)
    enrich_dfs = list()
    for enrich_file in enrich_dir.glob('*ALL.*.enrichment.txt'):
        enrich_df = pd.read_csv(enrich_file, sep='\t')
        enrich_df = enrich_df.loc[:, ENRICH_COL[enrich_type]]
        enrich_df.columns = ['term', 'qvalue']
        enrich_df.loc[:, 'compare'] = pattern.match(
            enrich_file.name).groups()[0]
        enrich_dfs.append(enrich_df)
    m_enrich_df = pd.concat(enrich_dfs)
    m_enrich_df.sort_values('qvalue', inplace=True)
    m_enrich_df = m_enrich_df.reset_index().drop('index', axis=1)
    top_enrich_dict = defaultdict(dict)
    for idx in m_enrich_df.index:
        term, qvalue, compare = m_enrich_df.loc[idx]
        top_enrich_dict[term][compare] = -np.log10(qvalue)
        if len(top_enrich_dict.keys()) == number:
            break
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f'top.{enrich_type}.enrich.json'
    with open(out_file, 'w') as out_inf:
        json.dump(top_enrich_dict, out_inf)


if __name__ == '__main__':
    fire.Fire(top_enriched_items)
