import fire
import delegator
import pandas as pd
from pathlib import PurePath
from pybedtools import BedTool


R_BIN = '/public/software/R/R-3.5.1/executable/bin/'
PLOT_R = PurePath(__file__).parent / 'lnc_cis_dis.R'


def lnc_neighbor_dis(lnc_nb, chr_size, chr_window):
    nb_file = PurePath(lnc_nb)
    outdir = nb_file.parent
    nb_df = pd.read_csv(nb_file)
    nb_df.loc[:, "chrom"] = [each.split(':')[0] for each in nb_df.Locus]
    nb_df.loc[:, 'end'] = [each.split(':')[1].split(
        '|')[0].split('-')[1] for each in nb_df.Locus]
    nb_df.loc[:, 'start'] = [int(each.split(':')[1].split(
        '|')[0].split('-')[0]) - 1 for each in nb_df.Locus]
    nb_df = nb_df.loc[:, ['chrom', 'start', 'end',
                          'PCC', 'lncRNA_gene', 'partnerRNA_gene']]
    nb_df.drop_duplicates(inplace=True)
    gene_bed_file = outdir / 'lnc.neighbor.bed'
    nb_df.to_csv(
        gene_bed_file, sep='\t',
        header=False, columns=[
            'chrom', 'start', 'end',
            'PCC'], index=False)
    gene_bed = BedTool(str(gene_bed_file))
    window_bed = BedTool(chr_window)
    intersect_obj = window_bed.intersect(gene_bed, wo=True)
    intersect_df = intersect_obj.to_dataframe()
    region_pcc = intersect_df.groupby(['chrom', 'start', 'end'])[
        'thickStart'].mean()
    region_count = intersect_df.groupby(['chrom', 'start', 'end']).size()
    region_count.name = 'gene_count'
    region_pcc.name = 'pcc'
    region_inf = pd.concat([region_pcc, region_count], axis=1)
    region_inf_df = pd.DataFrame(region_inf)
    region_inf_file = outdir / 'region.lnc.pcg.pcc.csv'
    region_inf_df.to_csv(region_inf_file)
    plt_cmd = (f'{R_BIN}/Rscript {PLOT_R} '
               f'--lnc_cis {region_inf_file} '
               f'--chrom_size {chr_size} '
               f'--out_dir {outdir}')
    print(plt_cmd)
    delegator.run(plt_cmd)


if __name__ == '__main__':
    fire.Fire(lnc_neighbor_dis)
