import fire
import gtfparse
import re
import sys
import pandas as pd
from pathlib import PurePath
from Bio import SeqIO, Seq


def gene_presentative_fa(gffread_fa, gtf=None,
                         gene_tr_map=None, transdecode=False):
    if gtf is not None:
        gtf_df = gtfparse.read_gtf(gtf)
        tr_gtf_df = gtf_df[
            gtf_df.feature == 'exon'].loc[
                :, ['transcript_id', 'gene_id']]
        tr_gtf_df = tr_gtf_df.drop_duplicates().set_index('transcript_id')
    elif gene_tr_map is not None:
        tr_gtf_df = pd.read_csv(gene_tr_map, header=None,
                                index_col=1, names=['gene_id'],
                                sep='\t')
    else:
        tr_gtf_df = None
    gene_fa_dict = dict()
    for seq_record in SeqIO.parse(gffread_fa, "fasta"):
        # replace dot in seq, diamond mkindex fail
        seq_record.seq = Seq.Seq(str(seq_record.seq).replace('.', '*'))
        tr_id = seq_record.id
        if transdecode:
            tr_id = re.match('(\S+).p\w+$', tr_id).groups()[0]
        if tr_gtf_df is None:
            if re.search('gene=(\S+)', seq_record.description):
                gene_id = re.search(
                    'gene=(\S+)', seq_record.description).groups()[0]
            else:
                sys.exit(
                    'Wrong fasta format. [>transcript_id gene=gene_id]: {}'.format(
                        seq_record.description
                    ))
        else:
            gene_id = tr_gtf_df.loc[tr_id].gene_id
        seq_len = len(seq_record.seq)
        seq_record.id = gene_id
        seq_record.description = ''
        if gene_id in gene_fa_dict:
            if gene_fa_dict[gene_id][0] >= seq_len:
                continue
        gene_fa_dict[gene_id] = [seq_len, seq_record]

    seq_record_list = []
    for gene_i in gene_fa_dict:
        seq_record_list.append(gene_fa_dict[gene_i][1])
    gffread_presentative_fa = PurePath(gffread_fa).with_suffix('.gene.fa')
    SeqIO.write(seq_record_list, gffread_presentative_fa, 'fasta')


if __name__ == '__main__':
    fire.Fire(gene_presentative_fa)
