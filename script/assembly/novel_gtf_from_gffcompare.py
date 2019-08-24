import fire
from HTSeq import GFF_Reader


NOVEL_TR_CODE = ('u', 'p')
OUT_ATTR = ('transcript_id', 'gene_id')


def novel_gtf(compare_gtf, outfile):
    novel_gtf_dict = dict()
    outfile_inf = open(outfile, 'w')
    for record in GFF_Reader(compare_gtf):
        if 'class_code' in record.attr:
            if record.attr['class_code'] in NOVEL_TR_CODE:
                novel_gtf_dict[record.attr['transcript_id']] = 1
                record.attr = {key: val for key,
                               val in record.attr.items()
                               if key in OUT_ATTR}
            else:
                continue
        elif record.attr['transcript_id'] in novel_gtf_dict:
            pass
        else:
            continue
        outline = record.get_gff_line().strip()
        outfile_inf.write(f'{outline};\n')
    outfile_inf.close()


if __name__ == '__main__':
    fire.Fire(novel_gtf)
