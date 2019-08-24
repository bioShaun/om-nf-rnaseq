import click
import pandas as pd
import os
import glob
from functools import reduce

CLICK_CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


class MutuallyExclusiveOption(click.Option):
    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(kwargs.pop('mutually_exclusive', []))
        help = kwargs.get('help', '')
        if self.mutually_exclusive:
            ex_str = ', '.join(self.mutually_exclusive)
            kwargs['help'] = help + (
                ' NOTE: This argument is mutually exclusive with '
                ' arguments: [' + ex_str + '].'
            )
        super(MutuallyExclusiveOption, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.mutually_exclusive.intersection(opts) and self.name in opts:
            raise click.UsageError(
                "Illegal usage: `{}` is mutually exclusive with "
                "arguments `{}`.".format(
                    self.name,
                    ', '.join(self.mutually_exclusive)
                )
            )

        return super(MutuallyExclusiveOption, self).handle_parse_result(
            ctx,
            opts,
            args
        )


def merge_files(df_list, df0=None, method='left'):
    if df0 is None:
        df0 = df_list.pop(0)
        return merge_files(df_list, df0, method=method)
    elif df_list:
        df1 = df_list.pop(0)
        new_df = pd.merge(df0, df1,
                          left_index=True,
                          right_index=True,
                          how=method)
        new_df.index.name = df0.index.name
        return merge_files(df_list, df0=new_df, method=method)
    else:
        return df0


@click.command(context_settings=CLICK_CONTEXT_SETTINGS)
@click.argument(
    'file_and_columns',
    nargs=-1)
@click.argument(
    'output',
    nargs=1)
@click.option(
    '-nh',
    '--noheader',
    is_flag=True,
    help='Input tables without header.',
)
@click.option(
    '-bn',
    '--by_colname',
    is_flag=True,
    help='Merge tables by identical columns.',
    cls=MutuallyExclusiveOption,
    mutually_exclusive=['by_row'],
)
@click.option(
    '-bh',
    '--by_row',
    is_flag=True,
    help='Merge tables by row.',
    cls=MutuallyExclusiveOption,
    mutually_exclusive=['by_colname'],
)
@click.option(
    '-na',
    '--na_rep',
    default='0',
    help='Value to replace NA in merged table. [default: 0]'
)
@click.option(
    '--method',
    default='left',
    type=click.Choice(['left', 'right', 'outer', 'inner'])
)
@click.option(
    '--out_sep',
    help='output file seperator, default is tab.',
    default="\t"
)
def main(file_and_columns, output, noheader, na_rep,
         by_colname, by_row, method, out_sep):
    '''Merge Multi files by matched column
    '''
    # seperate files and columns
    if by_colname or by_row:
        table_dfs = [pd.read_table(each_file)
                     for each_file in file_and_columns]
        if by_colname:
            merged_df = reduce(pd.merge, table_dfs)
        else:
            merged_df = pd.concat(table_dfs)
        merged_df.to_csv(output, sep=out_sep,
                         na_rep=na_rep, float_format='%.3f',
                         index=False)
        return 1

    files = list()
    columns = list()
    for n, each in enumerate(file_and_columns):
        if '*' in each:
            each_files = glob.glob(each)
            files.extend(each_files)
        elif os.path.isfile(each):
            files.append(each)
        else:
            break
    file_num = len(files)
    para_num = len(file_and_columns)
    if n < para_num - 1:
        columns = [int(each) for each in file_and_columns[n:]]
    # simple check
    if file_num < 2:
        click.secho('File number must more than 2.', color='red')
        return 1
    else:
        if columns:
            if file_num != len(columns):
                click.secho(
                    'File number should equal to column number.', color='red')
                return 1
        else:
            columns = [1 for each in range(file_num)]
        # read files based on merging columns
        # if tables has header
        if noheader:
            df_list = [
                pd.read_table(files[i], index_col=columns[i] - 1,
                              header=None)
                for i in range(file_num)
            ]
        else:
            df_list = [
                pd.read_table(files[i], index_col=columns[i] - 1)
                for i in range(file_num)
            ]
        # merge files
        merged_df = merge_files(df_list, method=method)
        header = not noheader
        merged_df.to_csv(output, sep=out_sep, header=header,
                         na_rep=na_rep, float_format='%.3f')
    return 1


if __name__ == '__main__':
    main()
