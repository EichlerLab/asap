#!/usr/bin/env python3
"""
Script to identify inheritance
Usage: ./process_DC_output.py --merge_map --disc_curve variants_bg_ALL-sv_insdel.tsv.gz -t 10
Author: Mei Wu https://github.com/projectoriented
"""
import pandas as pd
import numpy as np

import re
from itertools import groupby, repeat
import argparse
import sys

import multiprocessing as mp


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    parser.add_argument('--genotypes', metavar='gt_sv_insdel.tsv.gz', required=False,
                        help='Output from data table', default=None)
    parser.add_argument('--merge_map', metavar='merge_map_sv_insdel.tsv.gz', required=True,
                        help='Output from data table')
    parser.add_argument('--disc_curve', required=True, help='Output from discovery curve',
                        default=False)
    parser.add_argument('-t', required=False, type=int, help='Number of threads', default=1)
    parser.add_argument('--output', '-o', type=str, required=False, default=sys.stdout,
                        help='Output')
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Genotypes DF
    gt_df = pd.read_table(args.genotypes, index_col=['ID']) if args.genotypes else pd.DataFrame(index=["ID"])

    # Merge map DF
    m_df = pd.read_table(args.merge_map, index_col=['ID'])
    m_df.fillna('NONE', inplace=True)

    # Discovery curve DF
    dc_df = pd.read_table(args.disc_curve)

    # Sanity check desired output size
    print(f'Untouched discovery curve (rows, columns): {dc_df.shape}')

    # Subset to our desired samples so the working dataframe is smaller.
    dc_df = dc_df[~dc_df.MERGE_SAMPLES.str.contains('HG|NA')]

    n_chunks = dc_df.shape[0] / 2
    if args.t > 1:
        chunked_df = np.array_split(dc_df, n_chunks)
        pool = mp.Pool(args.t)
        with pool as p:
            result_df = p.starmap(combine_tables, zip(chunked_df, repeat(m_df), repeat(gt_df)))
            p.close()
            p.join()

        result_df = pd.concat(result_df)
        result_df.fillna('N/A', inplace=True)
    else:
        result_df = combine_tables(dc_df=dc_df, m_df=m_df, gt_df=gt_df)

    asd_parser_object = ASDParser(table_file_path=args.disc_curve)

    asd_parser_object.add_column('MERGE_VARIANTS', 'N/A')
    asd_parser_object.add_column('MERGE_GT', 'N/A')

    asd_parser_object.df.update(result_df)

    asd_parser_object.get_inheritance()

    # Sanity check desired output size
    print(f'Final output shape (rows, columns): {asd_parser_object.df.shape}')

    # Write out
    asd_parser_object.write_out(out_name=args.output)


class ASDParser:

    def __init__(self, table_file_path: str):
        self.df = pd.read_table(table_file_path).set_index('ID', drop=False)

    @staticmethod
    def _n_inherited(row: str) -> int:
        """
        Get the number of inherited events per family in a comma separated string.
        :param row: MERGE_SAMPLES column, e.g. '12456_s1,14455_fa,12456_fa,14455_p1'
        :return: An integer indicating number of events inherited.
        """
        starting_list = row.split(',')
        starting_list.sort()

        # Group the samples by family
        grouped = groupby(starting_list, lambda x: re.match(r'(.+)[-_]{1}(.+)', x).group(1))
        final_list = [','.join(g) for k, g in grouped]

        inherited = 0
        for str_entry in final_list:
            n_children = len(re.findall('s|p', str_entry))
            n_parents = len(re.findall('mo|fa', str_entry))
            if n_children >= 1 and n_parents >= 1:
                inherited += n_children
        return inherited

    @staticmethod
    def _inherited_by(row) -> str:
        """
        Get the inheritance for the event. Dataframe must contain 'FAMILY' and 'MERGE_SAMPLES' column.
        """
        member_conv = {
            'mo': 'MOTHER',
            'fa': 'FATHER'
        }

        # # Check if MERGE_GT is in the dataframe
        # if 'MERGE_GT' in row.index:
        #     # Get the index of the non-parent
        #     non_parent_idx = [i for i, x in enumerate(row.MERGE_SAMPLES.split(',')) if ("_p" in x) or ("_s" in x)][0]
        #
        #     # Use the index to get the genotype
        #     gt = row.MERGE_GT.split(',')[non_parent_idx]
        #
        #     if gt == '0|1':
        #         return 'MOTHER'
        #     elif gt == '1|0':
        #         return 'FATHER'
        #     elif gt == '1|1':
        #         return 'BOTH'
        #     else:
        #         return 'N/A'

        pattern = r'([mofa]{2})'
        fam_list = row.FAMILY.split(',')

        if (len(fam_list) == 1) and ('' not in fam_list):
            match = re.findall(pattern, row.MERGE_SAMPLES)
            if len(match) == 1:
                return member_conv.get(match[0])
            elif len(match) == 2:
                return 'BOTH'
            else:
                return 'N/A'
        else:
            return 'N/A'

    def get_asd_frequency(self):
        self.df['ASD_FREQ'] = self.df.MERGE_SAMPLES.apply(
            lambda x: len(re.findall(r'[_-]', x)))
        return self.df

    def get_pop_frequency(self):
        self.df['POP_FREQ'] = self.df.MERGE_SAMPLES.apply(
            lambda x: len(re.findall(r'NA|HG', x)))
        return self.df

    def get_families(self):
        self.df['FAMILY'] = self.df.apply(
            lambda x: ','.join(
                set(x[0] for x in re.findall(r'(\w+)[_-](\w+)', x.MERGE_SAMPLES))) if x.POP_FREQ == 0 else 'N/A',
            axis=1)
        return self.df

    def get_inherited_from(self):
        self.df['INHERITED_FROM'] = self.df.apply(
            lambda x: self._inherited_by(x) if x.POP_FREQ == 0 else 'N/A',
            axis=1
        )
        return self.df

    def stratify_children_counts(self):
        self.df = self.get_asd_frequency()
        self.df = self.get_pop_frequency()

        self.df['N_PROBAND'] = self.df.apply(
            lambda x: len(re.findall(r'[-_]p', x.MERGE_SAMPLES)) if (x.ASD_FREQ > 0) and (x.POP_FREQ == 0) else 0,
            axis=1)
        self.df['N_SIBLING'] = self.df.apply(
            lambda x: len(re.findall(r'[-_]s', x.MERGE_SAMPLES)) if (x.ASD_FREQ > 0) and (x.POP_FREQ == 0) else 0,
            axis=1)

        def _get_person(row):
            if row.POP_FREQ == 0:
                if row.N_PROBAND == 1 and row.N_SIBLING == 0:
                    return 'PROBAND'
                elif row.N_PROBAND == 0 and row.N_SIBLING == 1:
                    return 'SIBLING'
                else:
                    return 'N/A'
            else:
                return 'N/A'

        self.df['PERSON'] = self.df.apply(lambda x: _get_person(x), axis=1)
        return self.df

    def get_event_type(self):
        self.df = self.stratify_children_counts()

        self.df['PUTATIVE_DENOVO'] = self.df.apply(
            lambda x: 1 if (x.POP_FREQ == 0) and (x.ASD_FREQ == 1) and (sum([x.N_PROBAND, x.N_SIBLING]) == 1) else 0, axis=1)

        self.df['N_INHERITED'] = self.df.apply(
            lambda x: self._n_inherited(x.MERGE_SAMPLES) if (x.POP_FREQ == 0) else 0,
            axis=1)

        return self.df

    def get_inheritance(self):
        self.df = self.get_event_type()

        self.df = self.get_families()
        self.df = self.get_inherited_from()
        return self.df

    def add_column(self, col_name, val):
        self.df[col_name] = val
        return self.df

    def get_asd_subset(self):
        return self.df[(self.df.ASD_FREQ > 0) & (self.df.POP_FREQ == 0)]

    def write_out(self, out_name):
        self.df.to_csv(out_name, index=False, header=True, sep='\t')
        print(f'{out_name} written, bye.')


def id_in_sample(id_name: str, df):
    if id_name in df.index:
        df = df[df.index == id_name]
        samples = []
        for idx_label, col_entries in df.items():
            if col_entries[0] != 'NONE':
                samples.append((col_entries.name, col_entries[0]))

        sample_names = ','.join([x[0] for x in samples])
        variants = ','.join([x[1] for x in samples])

        del samples
    else:
        sample_names = ['N/A']
        variants = ['N/A']

    return pd.DataFrame(
        data={
            'MERGE_SAMPLES': [sample_names],
            'MERGE_VARIANTS': [variants]
        }, index=[id_name]
    )


def gt_in_sample(id_name: str, df):
    desired_genotypes = ['0|1', '1|0', '1|1']

    if id_name in df.index:
        df = df[df.index == id_name]
        samples = []
        for idx_label, col_entries in df.items():
            if col_entries[0] in desired_genotypes:
                samples.append((col_entries.name, col_entries[0]))

        genotypes = ','.join([x[1] for x in samples])

        del samples
    else:
        genotypes = ['N/A']

    return pd.DataFrame(
        data={
            'MERGE_GT': [genotypes]
        }, index=[id_name]
    )


# Combine the input tables
def combine_tables(dc_df, m_df, gt_df):
    result_gt_df = pd.DataFrame()
    result_id_df = pd.DataFrame()

    for idx, row in dc_df.iterrows():
        one_id = row.ID
        arr_gt_df = gt_in_sample(id_name=one_id, df=gt_df)
        arr_id_df = id_in_sample(id_name=one_id, df=m_df)
        result_id_df = pd.concat([result_id_df, arr_id_df])
        result_gt_df = pd.concat([result_gt_df, arr_gt_df])

    result_df = pd.concat([result_id_df, result_gt_df], axis=1)
    return result_df


if __name__ == "__main__":
    sys.exit(main())
