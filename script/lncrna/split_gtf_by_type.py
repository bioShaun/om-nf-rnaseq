import fire
import gtfparse
from pathlib import Path


GENCODE_CATEGORY_MAP = {
    'IG_C_gene': 'protein_coding',
    'IG_D_gene': 'protein_coding',
    'IG_J_gene': 'protein_coding',
    'IG_V_gene': 'protein_coding',
    'IG_LV_gene': 'protein_coding',
    'TR_C_gene': 'protein_coding',
    'TR_J_gene': 'protein_coding',
    'TR_V_gene': 'protein_coding',
    'TR_D_gene': 'protein_coding',
    'TEC': 'protein_coding',
    'nonsense_mediated_decay': 'protein_coding',
    'non_stop_decay': 'protein_coding',
    'retained_intron': 'lncRNA',
    'protein_coding': 'protein_coding',
    'ambiguous_orf': 'lncRNA',
    'Mt_rRNA': 'ncRNA',
    'Mt_tRNA': 'ncRNA',
    'miRNA': 'ncRNA',
    'misc_RNA': 'ncRNA',
    'rRNA': 'ncRNA',
    'snRNA': 'ncRNA',
    'snoRNA': 'ncRNA',
    'ribozyme': 'ncRNA',
    'sRNA': 'ncRNA',
    'scaRNA': 'ncRNA',
    'scRNA': 'ncRNA',
    'non_coding': 'lncRNA',
    'known_ncrna': 'ncRNA',
    '3prime_overlapping_ncrna': 'lncRNA',
    '3prime_overlapping_ncRNA': 'lncRNA',
    'vaultRNA': 'ncRNA',
    'processed_transcript': 'lncRNA',
    'lincRNA': 'lncRNA',
    'macro_lncRNA': 'lncRNA',
    'sense_intronic': 'lncRNA',
    'sense_overlapping': 'lncRNA',
    'antisense': 'lncRNA',
    'antisense_RNA': 'lncRNA',
    'bidirectional_promoter_lncRNA': 'lncRNA',
    'IG_pseudogene': 'pseudogene',
    'IG_D_pseudogene': 'pseudogene',
    'IG_C_pseudogene': 'pseudogene',
    'IG_J_pseudogene': 'pseudogene',
    'IG_V_pseudogene': 'pseudogene',
    'TR_V_pseudogene': 'pseudogene',
    'TR_J_pseudogene': 'pseudogene',
    'Mt_tRNA_pseudogene': 'pseudogene',
    'tRNA_pseudogene': 'pseudogene',
    'snoRNA_pseudogene': 'pseudogene',
    'snRNA_pseudogene': 'pseudogene',
    'scRNA_pseudogene': 'pseudogene',
    'rRNA_pseudogene': 'pseudogene',
    'misc_RNA_pseudogene': 'pseudogene',
    'miRNA_pseudogene': 'pseudogene',
    'pseudogene': 'pseudogene',
    'processed_pseudogene': 'pseudogene',
    'polymorphic_pseudogene': 'pseudogene',
    'retrotransposed': 'pseudogene',
    'transcribed_processed_pseudogene': 'pseudogene',
    'transcribed_unprocessed_pseudogene': 'pseudogene',
    'transcribed_unitary_pseudogene': 'pseudogene',
    'translated_processed_pseudogene': 'pseudogene',
    'translated_unprocessed_pseudogene': 'pseudogene',
    'unitary_pseudogene': 'pseudogene',
    'unprocessed_pseudogene': 'pseudogene',
    'novel_lncRNA': 'lncRNA',
    'TUCP': 'TUCP',
    'lncRNA': 'lncRNA'
}


def simplify_gene_type(gene_type):
    if gene_type in GENCODE_CATEGORY_MAP:
        sim_type = GENCODE_CATEGORY_MAP.get(gene_type)
        if sim_type == 'lncRNA':
            sim_type = f'annotated_{sim_type}'
        elif sim_type == 'ncRNA':
            sim_type = f'other_{sim_type}'
        else:
            pass
        return sim_type
    else:
        raise ValueError(gene_type)


def dfline2gtfline(dfline):
    basic_inf = dfline[:8]
    basic_inf.fillna('.', inplace=True)
    basic_inf.frame = '.'
    basic_inf_list = [str(each) for each in basic_inf]
    basic_inf_line = '\t'.join(basic_inf_list)
    attr_inf = dfline[8:]
    attr_inf_list = []
    for key, val in attr_inf.items():
        if val:
            attr_inf_list.append(f'{key} "{val}";')
    attr_inf_line = ' '.join(attr_inf_list)
    return f'{basic_inf_line}\t{attr_inf_line}\n'


def split_gtf(gtf, outdir, novel=False):
    gtf_df = gtfparse.read_gtf(gtf)
    if 'gene_type' in gtf_df.columns:
        gtf_df.loc[:, 'gene_biotype'] = gtf_df.gene_type
        gtf_df.drop('gene_type', axis=1, inplace=True)
    elif 'gene_biotype' in gtf_df.columns:
        pass
    else:
        gtf_df.loc[:, 'gene_biotype'] = 'protein_coding'

    type_label = 'gene_biotype'

    if novel:
        gtf_df.loc[
            :, type_label] = gtf_df.loc[:, type_label].map(
                GENCODE_CATEGORY_MAP)
    else:
        gtf_df.loc[
            :, type_label] = gtf_df.loc[:, type_label].map(
            simplify_gene_type)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for gt, grp in gtf_df.groupby(type_label):
        gt_file = outdir / f'{gt}.gtf'
        with open(gt_file, 'w') as gt_inf:
            for idx in grp.index:
                outline = dfline2gtfline(grp.loc[idx])
                gt_inf.write(outline)


if __name__ == '__main__':
    fire.Fire(split_gtf)
