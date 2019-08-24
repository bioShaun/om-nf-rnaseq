import jinja2
from bs4 import BeautifulSoup
import fire
from pathlib import PurePath, Path
import pandas as pd
import numpy as np
import delegator
from loguru import logger
import sys

BLAST_COLS = [
    'query',
    'identity',
    'length',
    'mismatch',
    'gapopen',
    'qstart',
    'qend',
    'sstart',
    'send',
    'evalue',
    'bitscore',
]

CURRENT_DIR = PurePath(__file__).parent
JS_FILE = CURRENT_DIR / 'jquery.imagemapster.min.js'
env = jinja2.Environment(loader=jinja2.FileSystemLoader(
    searchpath=str(CURRENT_DIR)))
template = env.get_template('pathway_template.html')
logger.add(sys.stderr,
           format="{time} {level} {message}",
           filter="my_module",
           level="INFO")


def fix_mouseon_over_tag(tag):
    if 'onmouseout' not in tag.attrs:
        tag['onmouseout'] = ""
    if 'onmouseover' not in tag.attrs:
        tag['onmouseover'] = ""
    else:
        tag['onmouseover'] = tag['onmouseover'].replace('"', "&quot;")
    return tag


def fix_url(tag):
    tag['href'] = f'https://www.genome.jp{tag["href"]}'
    return tag


def is_diff_gene(gene_id, diff_df, diff_cutoff):
    if gene_id not in diff_df.index:
        return False
    logfc, qvalue = diff_df.loc[gene_id, ['logFC', 'FDR']]
    cond1 = qvalue <= diff_cutoff['qvalue_cutoff']
    cond2 = np.abs(logfc) >= diff_cutoff['logfc_cutoff']
    if cond1 and cond2:
        if logfc > 0:
            return 'up'
        else:
            return 'down'
    else:
        return False


def new_title_stats(tag, blasttab_df, diff_df, diff_cutoff):
    genes = tag['href'].split('?')[-1].split('+')
    diff_stats = []
    titile_inf = tag['title'].split(", ")
    id_labes = set()
    for n, each_gene in enumerate(genes):
        each_analysis_gene = '--'
        diff_stat = '--'
        ori_gene = titile_inf[n]
        if each_gene in blasttab_df.index:
            each_analysis_gene = blasttab_df.loc[each_gene, 'query']
            diff_label = is_diff_gene(each_analysis_gene, diff_df, diff_cutoff)
            if diff_label:
                diff_stat = np.round(diff_df.loc[each_analysis_gene, "logFC"],
                                     decimals=3)
                id_labes.add(diff_label)
            else:
                diff_stat = 'Constant'
        diff_stats.append(f'{each_analysis_gene}|{ori_gene}|{diff_stat}')
    return diff_stats, id_labes


def title_id_attr(tags, blasttab_df, diff_df, diff_cutoff):
    new_tags = []
    for each_tag in tags:
        diff_stats, id_labes = new_title_stats(each_tag, blasttab_df, diff_df,
                                               diff_cutoff)
        if id_labes:
            if len(id_labes) == 1:
                each_tag['id'] = id_labes.pop()
            else:
                each_tag['id'] = 'mixed'
            each_tag['title'] = '\n'.join(diff_stats)
        else:
            each_tag['id'] = ""
        new_tags.append(each_tag)
    return new_tags


def displaty_area_tags(raw_html, blasttab_df, diff_df, diff_cutoff):
    with open(raw_html) as raw_html_inf:
        soup = BeautifulSoup(raw_html_inf, 'lxml')
    area_tags = soup.findAll('area')
    area_tags = map(fix_mouseon_over_tag, area_tags)
    area_tags = map(fix_url, area_tags)
    area_tags = title_id_attr(area_tags, blasttab_df, diff_df, diff_cutoff)
    return area_tags


def pathway_html(raw_html,
                 blasttab_df,
                 diff_df,
                 outfile,
                 logfc_cutoff=1,
                 qvalue_cutoff=0.05):
    diff_cutoff = {
        'logfc_cutoff': logfc_cutoff,
        'qvalue_cutoff': qvalue_cutoff
    }
    raw_html = Path(raw_html)
    if not raw_html.exists():
        logger.warning(f'{raw_html.name} not exists!')
        return False

    display_dictionary = {}
    # render pathway id
    display_dictionary['pathway_id'] = raw_html.stem
    # render area tags
    display_dictionary['area_tags'] = displaty_area_tags(
        raw_html, blasttab_df, diff_df, diff_cutoff)

    # output html file
    display_html = template.render(display_dictionary)
    with open(outfile, 'w') as out_inf:
        out_inf.write(display_html)

    return True


def pathway_results(pathway_db,
                    kegg_table,
                    blasttab,
                    diff_table,
                    outdir,
                    kegg_abbr,
                    logfc_cutoff=1,
                    qvalue_cutoff=0.05):
    kegg_df = pd.read_csv(kegg_table, sep='\t')
    diff_df = pd.read_csv(diff_table, index_col=0, sep='\t')
    blasttab_df = pd.read_csv(blasttab,
                              header=None,
                              index_col=1,
                              names=BLAST_COLS,
                              sep='\t')
    pathway_db = PurePath(pathway_db)
    outdir = Path(outdir)
    src_dir = outdir / 'src'
    src_dir.mkdir(parents=True, exist_ok=True)

    # generate html files
    raw_htmls = [
        pathway_db / f'html/{kegg_abbr}/{pathway}.html'
        for pathway in kegg_df.ID
    ]
    pngs = [
        pathway_db / f'png/{kegg_abbr}/{pathway}.png' for pathway in kegg_df.ID
    ]
    for n, raw_html in enumerate(raw_htmls):
        outfile = outdir / f'{raw_html.name}'
        if pathway_html(raw_html,
                        blasttab_df,
                        diff_df,
                        outfile,
                        logfc_cutoff=logfc_cutoff,
                        qvalue_cutoff=qvalue_cutoff):
            # cp pngs
            cp_cmd = f'cp {pngs[n]} {src_dir}'
            delegator.run(cp_cmd)
    delegator.run(f'cp {JS_FILE} {src_dir}')


if __name__ == '__main__':
    fire.Fire(pathway_results)
