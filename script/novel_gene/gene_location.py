import pandas as pd
import gtfparse
import fire


def gene_loc_inf(gtf_file, outfile):
    gtf_df = gtfparse.read_gtf(gtf_file)
    gene_chrom = gtf_df.groupby(['gene_id'])['seqname'].first()
    gene_start = gtf_df.groupby(['gene_id'])['start'].min()
    gene_end = gtf_df.groupby(['gene_id'])['end'].max()
    gene_strand = gtf_df.groupby(['gene_id'])['strand'].first()
    locus = pd.concat([gene_chrom, gene_start,
                       gene_end, gene_strand], axis=1)

    locus = pd.DataFrame(locus)
    locus.loc[~locus.strand.isin(['+', '-']), 'strand'] = '.'
    locus = locus.astype('str')
    locus.loc[:, 'Locus'] = locus.seqname.str.cat(locus.start, sep=':')
    locus.loc[:, 'Locus'] = locus.Locus.str.cat(locus.end, sep='-')
    locus.loc[:, 'Locus'] = locus.Locus.str.cat(locus.strand, sep='|')
    locus.to_csv(outfile, sep='\t', columns=['Locus'])


if __name__ == '__main__':
    fire.Fire(gene_loc_inf)
