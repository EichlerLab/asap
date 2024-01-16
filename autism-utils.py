#!/usr/bin/env python3
"""
Usage: ./autism-utils.py sample_info.tab --which_one sample_summary
"""
import argparse
import sys
from datetime import date

import pandas as pd

CURRENT_DATE = date.today().strftime('%Y%m%d')


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument('filepath', nargs=1)
    parser.add_argument('--which_one', choices=["sample_summary"])
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    filepath = args.filepath[0]
    which_one = args.which_one

    if which_one == "sample_summary":
        print(get_sample_summary(filepath=filepath))


def get_sample_summary(filepath):
    df = pd.read_table(filepath, header=0, names=['sample', 'sex', 'family_id', 'family_type', 'alt_id', 'family_role',
                                                  'cohort'])

    # Take away the parents.
    df = df[~df.family_role.str.contains("mo|fa")]

    def get_sex_pair(row):
        # Make sure [p]roband is always first and [s]ibling is second.
        sorted_row = row.sort_values(by=["sample"])
        sex_pair = "-".join(sorted_row.sex.tolist())
        return sex_pair

    df = df.groupby(["family_type", "family_id"]).apply(get_sex_pair).reset_index().rename(
        columns={0: "sex_pair"}).groupby(["sex_pair", "family_type"])["family_id"].count().reset_index()

    return df


if __name__ == "__main__":
    sys.exit(main())
