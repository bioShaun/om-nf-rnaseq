'''
get kegg pathview plots

'''

import sys
import os
import argparse
import envoy


PATHWAY_BASE = '/public/database/kegg_pathway/'

parser = argparse.ArgumentParser()
parser.add_argument(
    '--kegg_table',
    help='KOBAS kegg enrichment analysis result',
    required=True)
parser.add_argument(
    '--blast_out',
    help='ko blast result',
    required=True)
parser.add_argument(
    '--species',
    help='kegg species',
    required=True)
parser.add_argument(
    '--diff_out',
    help='diff analysis result',
    default='')
parser.add_argument(
    '--diff_method',
    help='diff analysis method',
    choices=['edgeR', 'DESeq2', 'Cuffdiff'],
    default='edgeR')
parser.add_argument(
    '--out_dir',
    help='output directory',
    required=True)
args = parser.parse_args()

script_dir = os.path.dirname(os.path.abspath(__file__))
pathview_script = os.path.join(script_dir, 'kegg_pathview_plot.R')
# get pathway database
kegg_pathway_dir = PATHWAY_BASE 
if not kegg_pathway_dir:
    kegg_pathway_dir = os.path.join(script_dir, 'kegg_pathway')
    sp_pathway_dir = os.path.join(kegg_pathway_dir, args.species)
    if not os.path.exists(sp_pathway_dir):
        os.makedirs(sp_pathway_dir)
        download_scirpt = os.path.join(script_dir, 'download_kegg.R')
        envoy.run('Rscript {download_scirpt} {args.species} {kegg_pathway_dir}'.format(
            **locals()
        ))
# exclude some problem pathway
exclude_pathway_id = ('01100', '01110')
problem_pathway_id = ('bta04215')


def invert_dict(d):
    return dict((v, k) for k, v in d.iteritems())


def get_fc_row(diff_method):
    if diff_method == 'edgeR':
        return 2
    elif diff_method == 'Cuffdiff':
        return 6
    elif diff_method == 'DESeq2':
        return 5
    else:
        sys.exit('wrong diff analysis method! : %s ' % diff_method)


def get_diff_fc_dict(diff_out, diff_method):
    fc_row = get_fc_row(diff_method)
    diff_fc_dict = {}
    with open(diff_out) as diff_out_info:
        for n, eachline in enumerate(diff_out_info):
            if n != 0:
                eachline_info = eachline.strip().split('\t')
                if eachline_info[fc_row - 1] != 'NA':
                    qurey_id = eachline_info[0]
                    fc = float(eachline_info[fc_row - 1])
                    diff_fc_dict[qurey_id] = fc
    return diff_fc_dict


def get_kegg_map(kegg_blast_out):
    kegg_map_dict = {}
    kegg_to_gene_dict = {}
    kegg_identity_dict = {}
    with open(kegg_blast_out) as kegg_blast_out_info:
        for eachline in kegg_blast_out_info:
            eachline_info = eachline.strip().split('\t')
            gene_id = eachline_info[0]
            kegg_gene = eachline_info[1].split(':')[1]
            identity = float(eachline_info[3])
            if kegg_gene not in kegg_to_gene_dict:
                kegg_to_gene_dict[kegg_gene] = gene_id
                kegg_identity_dict[kegg_gene] = identity
            else:
                if identity > kegg_identity_dict[kegg_gene]:
                    kegg_to_gene_dict[kegg_gene] = gene_id
                    kegg_identity_dict[kegg_gene] = identity
    kegg_map_dict = invert_dict(kegg_to_gene_dict)
    return kegg_map_dict


def plot_pathview(species, pathview_id, each_pathway_kegg_fc_out, database, out_dir):
    cmd = 'Rscript %s %s %s %s %s %s' % (
        pathview_script, species, pathview_id, each_pathway_kegg_fc_out, database, out_dir)
    envoy.run(cmd)
    os.system('rm %s' % each_pathway_kegg_fc_out)


def kegg_pathway_plot(kegg_out, kegg_map_dict, diff_fc_dict=dict()):
    with open(kegg_out) as kegg_out_info:
        for eachline in kegg_out_info:
            eachline_info = eachline.strip().split('\t')
            if (not eachline_info[0].startswith('#')) and len(eachline_info) == 9:
                pathway_id = eachline_info[2]
                pathview_id = pathway_id.split(args.species)[1]
                kegg_pic1 = os.path.join(
                    args.out_dir, '%s.pathview.png' % pathway_id)
                kegg_pic2 = os.path.join(
                    args.out_dir, '%s.pathview.gene.png' % pathway_id)
                if os.path.isfile(kegg_pic1) and os.stat(kegg_pic1).st_size and os.path.isfile(kegg_pic2) and os.stat(kegg_pic2).st_size:
                    pass
                elif pathview_id in exclude_pathway_id or pathway_id in problem_pathway_id:
                    target_pic = os.path.join(
                        kegg_pathway_dir, '%s/%s.png' % (args.species, pathway_id))
                    os.system('cp %s %s' % (target_pic, kegg_pic1))
                    os.system('cp %s %s' % (target_pic, kegg_pic2))
                else:
                    each_pathway_kegg_fc = os.path.join(
                        args.out_dir, '%s.keggid' % pathway_id)
                    each_pathway_kegg_fc_out = open(each_pathway_kegg_fc, 'w')
                    pathway_gene_list = eachline_info[7].split('|')
                    for each_gene in pathway_gene_list:
                        if each_gene in kegg_map_dict:
                            kegg_id = kegg_map_dict[each_gene]
                            if diff_fc_dict:
                                fc = diff_fc_dict[each_gene]
                            else:
                                fc = 1
                            each_pathway_kegg_fc_out.write(
                                '%s\t%s\n' % (kegg_id, fc))
                    each_pathway_kegg_fc_out.close()
                    plot_pathview(args.species, pathview_id,
                                  each_pathway_kegg_fc, kegg_pathway_dir, args.out_dir)


if __name__ == '__main__':
    diff_fc_dict = dict()
    if args.diff_out:
        diff_fc_dict = get_diff_fc_dict(args.diff_out, args.diff_method)
    kegg_map_dict = get_kegg_map(args.blast_out)
    kegg_pathway_plot(args.kegg_table, kegg_map_dict, diff_fc_dict)
