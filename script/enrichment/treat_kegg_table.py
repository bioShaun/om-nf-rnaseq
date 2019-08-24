import click
import os


@click.command(
    context_settings=dict(help_option_names=['-h', '--help'])
)
@click.option('-k', '--kegg',
              type=click.Path(exists=True),
              help='kegg enrich analysis output',
              required=True)
def treat_KEGG_table(kegg):
    kegg_out_dir, kegg_out_name = os.path.split(kegg)
    kegg_tmp_file = os.path.join(kegg_out_dir, 'tmp.%s' % kegg_out_name)
    os.system('mv %s %s' % (kegg, kegg_tmp_file))

    def check_KOBAS_out(kobas_out):
        kobas_out_info = open(kobas_out, 'r').readlines()
        flag = True
        for eachline in kobas_out_info:
            if not eachline.startswith("#") and len(eachline.strip().split('\t')) == 9:
                flag = True
                break
        else:
            flag = False
        return flag

    if check_KOBAS_out(kegg_tmp_file):
        kegg_out_info = open(kegg, 'w')
        with open(kegg_tmp_file, 'r') as kegg_tmp_file_info:
            count = 0
            for eachline in kegg_tmp_file_info:
                if len(eachline.strip().split('\t')) == 9:
                    if count == 0 and eachline.startswith("#"):
                        kegg_out_info.write(eachline)
                        count += 1
                    elif not eachline.startswith("#"):
                        kegg_out_info.write(eachline)
        kegg_out_info.close()
    os.system('rm %s' % (kegg_tmp_file))


if __name__ == '__main__':
    treat_KEGG_table()
