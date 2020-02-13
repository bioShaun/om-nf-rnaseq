import pandas as pd
import delegator
import fire
from pathlib import Path
from shutil import which
import sys
from loguru import logger
import re

logger.add(sys.stderr, format="{time} {level} {message}",
           filter="my_module", level="INFO")
# logger.add("file_{time}.log")

BASE_DIR = Path(__file__).parent

BLAST_HEADER = [
    'gene_id',
    'protein_id',
    'blast_identity',
    'blast_evalue',
    'bitscore',
    'stitle'
]

DB_LOCATION = {
    'swissprot': '/public/database/swissprot/uniprot_sprot.fasta'
}

TRANSDECODE = '/usr/bin/perl /public/software/TransDecoder/TransDecoder-v5.3.0/TransDecoder.LongOrfs'
GENE_PEP_PY = BASE_DIR / 'gene_presentative_fa.py'


def check_exe(program):
    '''check required progrom is in PATH'''
    if which(program) is not None:
        return 1
    else:
        logger.error('{program} is required for analysis.',
                     program=program)
        sys.exit(1)


def check_output(output, cmd, std_err=None):
    output = Path(output)
    if output.exists() and output.stat().st_size:
        logger.success(cmd)
        return True
    else:
        if std_err is not None:
            logger.error('{cmd}', cmd=cmd)
            logger.error(std_err)


def exe_cmd(program, cmd, output):
    '''execute shell cmd'''
    if check_output(output, cmd):
        return True
    else:
        check_exe(program)
        r = delegator.run(cmd)
        check_output(output, cmd, r.err)


def orfid2geneid(orfid, tr2gene_df):
    tr_id = '.'.join(orfid.split('.')[:-1])
    return tr2gene_df.loc[tr_id].gene_id


def extract_swissprot_inf(swissprot_des, stat):
    basic_stat = ('sp', 'des')
    if stat == 'des':
        re_pattern = re.compile('sp\|\w+\|\w+ (.*) OS=.*')
    elif stat == 'sp':
        re_pattern = re.compile('.* OS=(.*?) \w+=.*')
    elif stat == 'gn':
        re_pattern = re.compile('.* GN=(.*?) \w+=.*')
    else:
        logger.error('unsupported stat: {stat}',
                     stat=stat)
        sys.exit(1)
    if re_pattern.match(swissprot_des):
        stat_inf = re_pattern.match(swissprot_des).groups()[0]
        if stat in basic_stat and stat_inf is None:
            logger.error(
                'failed to find [{stat}] in swissprot description:{des}',
                stat=stat, des=swissprot_des)
            sys.exit(1)
        return stat_inf
    else:
        if stat in basic_stat:
            logger.error(
                'regrex pattern failed to match swissprot description:{des}',
                des=swissprot_des)
            sys.exit(1)
        else:
            return None


def swissprot_annotation(input_file, gene2tr,
                         input_type='fasta',
                         db='swissprot',
                         thread=40):
    '''annotate novel gene using swissprot database'''
    input_file = Path(input_file)
    # step1 orf annotation
    orf_cmd = f'{TRANSDECODE} -t {input_file} --gene_trans_map {gene2tr}'
    orf_file = f'{input_file}.transdecoder_dir/longest_orfs.pep'
    orf_cmd_err = exe_cmd('TransDecoder.LongOrfs', orf_cmd, orf_file)

    # step2 find longest orf for gene
    long_pep_cmd = f'python {GENE_PEP_PY} {orf_file} --gene_tr_map {gene2tr} --transdecode'
    long_pep_file = f'{input_file}.transdecoder_dir/longest_orfs.gene.fa'
    long_pep_cmd_err = exe_cmd('python', long_pep_cmd, long_pep_file)

    # step2 blast to swissprot database
    blast_out = input_file.with_suffix('.blasttab')
    blast_cmd = f'blastp -query {long_pep_file} -db {DB_LOCATION[db]} \
-out {blast_out} -evalue 1e-5 -max_target_seqs 1 -num_threads {thread} \
-outfmt "6 qseqid sseqid pident evalue bitscore stitle"'
    blast_cmd_err = exe_cmd('blastp', blast_cmd, blast_out)

    # step3 extract annotation from blast outfile
    blast_df = pd.read_table(blast_out, header=None,
                             names=BLAST_HEADER)
    best_match = blast_df.groupby(['gene_id'])['bitscore'].idxmax()
    best_blast_df = blast_df.reindex(best_match)
    best_blast_df.loc[:, 'species'] = [extract_swissprot_inf(each, 'sp')
                                       for each in best_blast_df.stitle]
    best_blast_df.loc[:, 'gene_name'] = [extract_swissprot_inf(each, 'gn')
                                         for each in best_blast_df.stitle]
    best_blast_df.loc[:, 'description'] = [extract_swissprot_inf(each, 'des')
                                           for each in best_blast_df.stitle]
    anno_file = input_file.with_suffix('.annotation.txt')
    best_blast_df.to_csv(anno_file, sep='\t', index=False,
                         columns=['gene_id', 'gene_name', 'gene_description'],
                         na_rep='--')
    anno_detail_file = input_file.with_suffix('.annotation.detail.txt')
    best_blast_df.to_csv(anno_detail_file, index=False,
                         columns=['gene_id', 'protein_id', 'blast_identity',
                                  'blast_evalue', 'species', 'gene_name',
                                  'gene_description'],
                         na_rep='--', sep='\t')


if __name__ == '__main__':
    fire.Fire(swissprot_annotation)
