#!/usr/bin/env python3

"""
Usage: ./check_denovo.py -s 14455_p1 -f size50_2_4.tsv.gz size50_2_4.tsv.gz size50_2_4.tsv.gz
"""

import argparse
import sys
import pandas as pd

TARGET_COL_NAMES = ['ORIGIN', 'DENOVO_CORRECTION', 'SUPPORT_COUNT', 'INHERITED_FROM_CORRECTION']


def check_inheritance(row, family_dict: dict):

    father = family_dict['father']
    mother = family_dict['mother']
    other_child = family_dict['other_child']

    n_parents = 0
    # n_child = 0
    one_parent = ''
    # if row[father] > 0 and row[mother] > 0 and row[other_child] == 0:
    #     n_parents = 2
    # elif row[father] > 0 and row[mother] == 0 and row[other_child] == 0:
    #     n_parents = 1
    #     one_parent = father
    # elif row[father] == 0 and row[mother] > 0 and row[other_child] == 0:
    #     n_parents = 1
    #     one_parent = mother

    if row[father] > 0 and row[mother] > 0:
        n_parents = 2
    elif row[father] > 0 and row[mother] == 0:
        n_parents = 1
        one_parent = father
    elif row[father] == 0 and row[mother] > 0:
        n_parents = 1
        one_parent = mother

    # if row[other_child] > 0:
    #     n_child = 1

    # # Get who it was inherited from
    # if n_child == 0:
    #     if n_parents == 2:
    #         return 'BOTH'
    #     elif n_parents == 1:
    #         return "MOTHER" if ("_mo" in one_parent) else "FATHER"
    #     else:
    #         return 'NONE'
    # else:
    #     return 'NONE'

    # Get who it was inherited from

    if n_parents == 2:
        return 'BOTH'
    elif n_parents == 1:
        return "MOTHER" if ("_mo" in one_parent) else "FATHER"
    else:
        return "NONE"


def verify_denovo(items: list, origin):
    sample = origin

    df = pd.concat([pd.read_csv(item,
                                sep='\t', header=0, index_col=['ID', 'SAMPLE']) for item in items])

    pt = pd.pivot_table(df, values=['SUPPORT_COUNT'], index=['ID'], columns=['SAMPLE'], fill_value=0)

    # Flatten the table to just one level
    pt.columns = pt.columns.get_level_values(1)

    # Add DENOVO information
    pt[TARGET_COL_NAMES[1]] = pt.loc[:, (~pt.columns.str.contains(sample))].apply(
        lambda x: 'YES' if x.sum() == 0 else 'NO', axis=1)

    # Get support count string
    pt[TARGET_COL_NAMES[2]] = pt.loc[:, ~pt.columns.str.contains(TARGET_COL_NAMES[1])].apply(
        lambda x: ';'.join(x.index + '=' + x.astype(str)), axis=1)

    family_dict = {
        "mother": pt.columns[pt.columns.str.contains("_mo")].tolist()[0],
        "father": pt.columns[pt.columns.str.contains("_fa")].tolist()[0],
        "other_child": pt.columns[~pt.columns.str.contains(f"_mo|_fa|{sample}", regex=True)].tolist()[0]
    }

    # Remove column name
    pt.columns.name = ''

    # Get who the event was inherited from
    pt[TARGET_COL_NAMES[3]] = pt.apply(lambda x: check_inheritance(row=x, family_dict=family_dict), axis=1)

    # ORIGIN
    pt[TARGET_COL_NAMES[0]] = sample

    return pt


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', metavar='14455_p1',
                        help='query/origin sample name')
    parser.add_argument('-f', metavar='p1.tsv.gz s1.tsv.gz fa.tsv.gz', nargs='+')
    parser.add_argument('--output', '-o', type=str, required=False, default=sys.stdout,
                        help='Output')
    return parser


def main():
    """Get bed"""
    parser = get_parser()
    args = parser.parse_args()

    pivot_table = verify_denovo(items=args.f, origin=args.s)

    inferred_sample_names = pivot_table.columns[pivot_table.columns.str.contains('[0-9]', regex=True)].to_list()
    out_file_header_order = TARGET_COL_NAMES + inferred_sample_names

    pivot_table.to_csv(args.output, index=True,
                       columns=out_file_header_order, sep='\t')


if __name__ == "__main__":
    sys.exit(main())
