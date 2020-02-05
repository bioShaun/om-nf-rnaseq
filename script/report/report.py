import fire
import glob
import yaml
import jinja2
import pandas as pd
from pathlib import Path, PurePath
from distutils import dir_util, file_util
from pybtex.database import parse_file
from pylatexenc.latex2text import LatexNodes2Text


# set global values
SCRIPT_DIR = Path(__file__).absolute().parent
TEMP_DIR = SCRIPT_DIR / 'template'
TEMP_IMG_DIR = 'img/analysis/'
CFG_FILE = SCRIPT_DIR / 'process.yml'
BIBTEX = SCRIPT_DIR / 'ref.bib'
REPORT_CFG = yaml.load(open(CFG_FILE), Loader=yaml.FullLoader)
BIB_DATA = parse_file(BIBTEX)

STAR_TABLE_COLS = [
    'Uniquely mapped reads %',
    'Mismatch rate per base, %',
    '% of reads mapped to multiple loci',
    '% of reads mapped to too many loci',
    '% of reads unmapped: too many mismatches',
    '% of reads unmapped: too short',
    '% of reads unmapped: other',
]

# jinja2 load template
env = jinja2.Environment(
    loader=jinja2.FileSystemLoader(
        searchpath=f'{TEMP_DIR}'
    )
)
template = env.get_template('index.html')


def format_reads_df(reads_df):
    int_df = reads_df.iloc[:, 0:4]
    float_df = reads_df.iloc[:, -4:]
    float_df = float_df.applymap(lambda x: f'{x:.3f}')
    int_df = int_df.astype('int').applymap(lambda x: f'{x:,}')
    clean_df = pd.concat([int_df, float_df], axis=1)
    clean_df.index.name = 'Item'
    return clean_df


def format_star_df(star_df):
    clean_df = star_df.loc[:, STAR_TABLE_COLS]
    clean_df.index.name = 'Item'
    return clean_df


TABLE_CLEANER = {
    'data_summary': format_reads_df,
    'mapping_summary': format_star_df,
}


def table2dict(table_file, sep='\t', format_func_map=TABLE_CLEANER):
    name = PurePath(table_file).stem
    format_func = format_func_map.get(name)
    table_dict = dict()
    if table_file.is_file():
        table_df = pd.read_csv(table_file, sep=sep, index_col=0)
        if format_func is not None:
            table_df = format_func(table_df)
        table_df.sort_index(inplace=True)
        table_df = table_df.reset_index()
        for idx_i in table_df.index:
            table_dict.setdefault(
                f'{name}_body', []).append(list(table_df.loc[idx_i]))
        table_dict[f'{name}_header'] = list(table_df.columns)
    return table_dict


def plotitem2report(result_dir, report_dir, plot_flag):
    plot_dict = dict()
    plot_item = result_dir / f'{plot_flag}'
    if plot_item.is_dir():
        outpath = report_dir / f'{TEMP_IMG_DIR}/{plot_flag}'
        outpath.mkdir(exist_ok=True, parents=True)
        plot_list = sorted(glob.glob(f'{plot_item}/*png'))
        if plot_list:
            for plot in plot_list:
                file_util.copy_file(plot, outpath)
            plot_dict[plot_flag] = [
                f'{TEMP_IMG_DIR}/{plot_flag}/{PurePath(each).name}'
                for each in plot_list]
    else:
        plot_name = plot_flag
        plot_flag = PurePath(plot_flag).stem
        plot_name = PurePath(plot_name).name
        outpath = report_dir / f'{TEMP_IMG_DIR}'
        if plot_item.is_file():
            file_util.copy_file(plot_item, outpath)
            plot_dict[plot_flag] = f'{TEMP_IMG_DIR}/{plot_name}'
    return plot_dict


def cite_item(item, sep):
    if item:
        return f'{item}{sep}'
    else:
        return ''


def cite_detail(cite_name):
    if cite_name not in BIB_DATA.entries:
        return None
    cite_obj = BIB_DATA.entries[cite_name]
    author_list = [str(each) for each in cite_obj.persons['author']]
    author = ', '.join(author_list) + '. '
    title = cite_obj.fields['title'] + '. '
    journal = cite_item(cite_obj.fields.get('journal'), '. ')
    volume = cite_item(cite_obj.fields.get('volume'), ':')
    page = cite_item(cite_obj.fields.get('pages'), ', ')
    if volume and page:
        index = f'{volume}{page}'
    else:
        index = ''
    year = cite_item(cite_obj.fields.get('year'), '.')
    cite_str = f'{author}{title}{journal}{index}{year}'
    cite_str = LatexNodes2Text().latex_to_text(cite_str)
    return cite_str


class ReportGenerator:

    def __init__(self, result_dir, report_dir):
        self.cite = 1
        self.p_order = REPORT_CFG['process_order']
        self.rs_dir = Path(result_dir)
        self.ro_dir = Path(report_dir)
        self._report_dict = {'cite': [], 'software': []}
        self.process_cite = None
        self.software_file = self.rs_dir / REPORT_CFG['software_file']
        self.software_name = REPORT_CFG['software_name']
        self.software_params = REPORT_CFG['software_params']

    def add_cite_software(self):
        software_df = pd.read_csv(self.software_file, index_col=0, sep='\t')
        for cite_label in self.process_cite:
            if cite_label in self._report_dict:
                for cite in self.process_cite[cite_label]:
                    cite_detail_str = cite_detail(cite)
                    if cite_detail_str:
                        self._report_dict['cite'].append(cite_detail_str)
                        self.cite += 1
                    software_name = self.software_name.get(cite, cite)
                    if software_name in software_df.index:
                        software_version = software_df.loc[
                            software_name].version
                        software_params = self.software_params.get(
                            cite, 'default')
                        self._report_dict['software'].append(
                            [software_name, software_version, software_params])

    @property
    def report(self):
        for process in self.p_order:
            process_cfg = REPORT_CFG['process'][process]
            process_label = self.rs_dir / process_cfg['label']
            process_table = process_cfg.get('table', [])
            process_plot = process_cfg.get('plot', [])
            self.process_cite = process_cfg.get('cite', [])
            if process_label.exists():
                self._report_dict[process] = True
                for table in process_table:
                    table_file = self.rs_dir / table
                    self._report_dict.update(
                        table2dict(table_file, sep=','))
                for plot in process_plot:
                    print(plot)
                    self._report_dict.update(
                        plotitem2report(self.rs_dir, self.ro_dir, plot)
                    )
            self._report_dict[f'{process}_cite'] = self.cite
            self.add_cite_software()
        return self._report_dict


def rnaseq_report(result_dir, proj_name='test', report_dir=None):
    result_dir = Path(result_dir)
    if report_dir is None:
        report_dir = result_dir / 'report'
    else:
        report_dir = Path(report_dir)
    if report_dir.is_dir():
        dir_util.remove_tree(report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    dir_util.copy_tree(f'{SCRIPT_DIR}/template/', f'{report_dir}')

    display_dictionary = {}
    display_dictionary['project_name'] = proj_name

    report_content = ReportGenerator(result_dir, report_dir)
    display_dictionary.update(report_content.report)

    display_html = template.render(display_dictionary)
    report_html = report_dir / 'index.html'
    with open(report_html, 'w', encoding="utf-8") as out_inf:
        out_inf.write(display_html)


if __name__ == "__main__":
    fire.Fire(rnaseq_report)
