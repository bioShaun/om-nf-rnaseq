import fire
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import PurePath


def get_cor(exp_df, gene1, gene2):
    try:
        return stats.pearsonr(exp_df.loc[gene1], exp_df.loc[gene2])
    except KeyError:
        return (np.nan, np.nan)


def lnc_pcg_neighbour_cor(gene_tpm_file, gene_neighbour_file,
                          grp_inf, pcg_feature):
    gene_tpm_df = pd.read_table(gene_tpm_file, index_col=0)
    grp_df = pd.read_csv(grp_inf, header=None, names=[
                         'sample_id', 'grp_id'], sep='\t')
    gf_df = pd.read_csv(pcg_feature)
    gf_df.drop(['Isoform_number', 'Length', 'Exon_number',
                'gene_biotype'], inplace=True, axis=1)
    gene_tpm_df = gene_tpm_df.loc[:, grp_df.sample_id]
    norm_gene_tpm_df = np.log2(gene_tpm_df + 1)
    gene_neighbour_df = pd.read_table(gene_neighbour_file)
    gene_neighbour_df = gene_neighbour_df.merge(
        gf_df, left_on='partnerRNA_gene',
        right_on='gene_id', how='left')
    gene_neighbour_df.drop(['gene_id'], inplace=True, axis=1)
    pearson_cor_list = []
    pearson_cor_p_list = []
    for row_i in gene_neighbour_df.index:
        lnc_i = gene_neighbour_df.loc[row_i, 'lncRNA_gene']
        pcg_i = gene_neighbour_df.loc[row_i, 'partnerRNA_gene']
        if pcg_i == '--':
            pearson_cor_list.append(np.nan)
            pearson_cor_p_list.append(np.nan)
            continue
        else:
            pcg_lnc_cor, pcg_lnc_cor_pval = get_cor(
                norm_gene_tpm_df, pcg_i, lnc_i)
            pearson_cor_list.append(pcg_lnc_cor)
            pearson_cor_p_list.append(pcg_lnc_cor_pval)
    gene_neighbour_df.loc[:, 'PCC'] = pearson_cor_list
    gene_neighbour_df.loc[:, 'PCC_Pvalue'] = pearson_cor_p_list
    gene_neighbour_df.dropna(inplace=True)
    gene_neighbour_file = PurePath(gene_neighbour_file)
    outfile_name = gene_neighbour_file.stem
    new_name = f'{outfile_name}_correlation.csv'
    outfile = PurePath(gene_neighbour_file).with_name(new_name)
    gene_neighbour_df.to_csv(outfile, index=False)
    best_nb = gene_neighbour_df[gene_neighbour_df.isBest == 1]
    best_dis = best_nb.groupby(['lncRNA_gene', 'partnerRNA_gene'])[
        'distance'].min()
    best_dis = pd.DataFrame(best_dis).reset_index()
    best_nb = best_nb.loc[:, ['lncRNA_gene', 'partnerRNA_gene', 'gene_name',
                              'gene_description', 'PCC', 'PCC_Pvalue']
                          ].drop_duplicates()
    best_nb = best_dis.merge(best_nb)
    best_nb_name = f'best_{outfile_name}_correlation.csv'
    best_nb_file = PurePath(gene_neighbour_file).with_name(best_nb_name)
    best_nb.to_csv(best_nb_file, index=False)


if __name__ == '__main__':
    fire.Fire(lnc_pcg_neighbour_cor)
