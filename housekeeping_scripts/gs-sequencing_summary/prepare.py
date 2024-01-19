#!/usr/bin/env python3
"""
Usage: ./prepare.py --tech nanopore --cohort clinical --proj_dir /path/to/long_read_archive/nobackups/
Author: Mei Wu, https://github.com/projectoriented
"""
import os
import glob
import sys
import logging
import argparse

import pandas as pd
import pysam

LOG = logging.getLogger()
logging.basicConfig(stream=sys.stderr, level="INFO", format='%(asctime)s - %(levelname)s - %(message)s')


def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    parser.add_argument('--outpath', required=False, type=str, help='Sheet name in the Google Sheet', default=sys.stdout)
    parser.add_argument('--tech', required=True, type=str, choices=["nanopore", "PacBio_HiFi"])
    parser.add_argument('--cohort', required=True, type=str, help='Sheet name in the Google Sheet')
    parser.add_argument('--proj_dir', required=True, type=str, help='Sheet name in the Google Sheet')

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    proj_dir = args.proj_dir
    cohort = args.cohort
    tech = args.tech
    outpath = args.outpath
    combine_quick_stats(prefix=proj_dir, cohort=cohort, tech=tech, outpath=outpath)


def get_sam_tag(fp, target_tag="SM") -> str:
    """
    Only supporting HiFi data right now.
    """

    if not os.path.exists(fp):
        return ""

    # Suppress I/O message on bam.bai not being there.
    save = pysam.set_verbosity(0)

    try:
        header_dict = pysam.AlignmentFile(fp, mode="rb", check_sq=False).header.to_dict()
    except OSError:
        header_dict = pysam.AlignmentFile(fp, mode="rb", check_sq=False, ignore_truncation=True).header.to_dict()
        LOG.info(f"{fp} is a truncated file.")

    pysam.set_verbosity(save)

    target_info = []
    for x in header_dict.keys():
        for items in header_dict.get(x):
            if isinstance(items, dict):
                value = items.get(target_tag, "")
                if value:
                    # Remove trailing bits.
                    target_info.append(value.strip())

    if target_info:
        if len(target_info) > 1:
            unique_list = list(set(target_info))
            return ",".join(unique_list)
        else:
            return target_info[0]
    else:
        return ""


def combine_quick_stats(prefix, cohort, outpath, tech="nanopore"):
    prefix = os.path.join(prefix, cohort, "*", "raw_data", tech)

    if tech == "nanopore":
        file_pattern = f'{prefix}/*/quick_stats/*.txt'
        target_columns = ["CELL", "LIBRARY", "BASECALLER", "BASECALLER_VERSION", "BASECALLER_MODEL", "FILE_PATH"]
    else:
        file_pattern = f'{prefix}/quick_stats/*.tsv'
        target_columns = ["CELL", "FILE_PATH"]

    quick_stats = glob.glob(file_pattern)
    if len(quick_stats) == 0:
        LOG.error(f"No files found for {file_pattern}")
        sys.exit()

    df = pd.concat(
        [
            pd.read_table(file, header=0, skipfooter=1, engine='python', dtype=str) for file in quick_stats
        ]
    ).rename(columns={"FILES": "FILE_PATH"})

    # Make sure we follow the real path and not soft links.
    df["is_link"] = df.apply(lambda row: 1 if os.path.islink(row.FILE_PATH) else 0, axis=1)
    df = df.query("is_link == 0").reset_index(drop=True).drop(columns=["is_link"])

    # Make sure we have the real path from hard links.
    for row in df.itertuples():
        try:
            if os.stat(row.FILE_PATH).st_nlink > 1:
                df.loc[row.Index, "true_path"] = os.path.realpath(row.FILE_PATH)
            else:
                df.loc[row.Index, "true_path"] = row.FILE_PATH
        except FileNotFoundError:
            df.drop([row.Index], inplace=True)
            LOG.info(f"{row.FILE_PATH} does not exist, skipping.")

    df = df.drop_duplicates(subset=["true_path"], keep="first").drop(columns=["true_path"])

    def add_new_columns(row):
        splitted_file_path = row.FILE_PATH.split(os.sep)

        sample = splitted_file_path[8]
        out_list = [sample]
        if tech == "nanopore":
            lib = splitted_file_path[11]
            basecaller = splitted_file_path[14]
            basecaller_version = splitted_file_path[15]
            basecaller_model = splitted_file_path[16]
            out_list.extend([lib, basecaller, basecaller_version, basecaller_model])
        return pd.Series(out_list)

    if tech == "nanopore":
        df[["SAMPLE", "LIBRARY", "BASECALLER", "BASECALLER_VERSION", "BASECALLER_MODEL"]] = df.apply(
            lambda row: add_new_columns(row=row), axis=1
        )
    else:
        # Add the sample name to the resulting file
        df["SAMPLE"] = df.apply(lambda row: add_new_columns(row=row), axis=1)

        # Add the target sam tag to the resulting file
        df["bam_fp"] = df.FILE_PATH.str.replace(r".fastq.gz$", ".bam", regex=True)
        df["RUN_SOURCE"] = df.apply(lambda row: get_sam_tag(fp=row.bam_fp, target_tag="SM"), axis=1)
        df.drop(columns=["bam_fp"], inplace=True)

    df.drop_duplicates(subset=target_columns, keep='last', inplace=True)
    df.reset_index(inplace=True, drop=True)

    df.to_csv(outpath, sep="\t", header=True, index=False)


if __name__ == "__main__":
    sys.exit(main())
