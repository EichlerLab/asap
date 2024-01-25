#!/usr/bin/env python3
"""
Usage: ./sex-estimator.py alignment.bam --threads 10 --prefix 11611_fa-GRCh38-HiFi-aln
Check the sex of sample.
Author: Mei Wu, https://github.com/projectoriented
"""
import os
import sys
import logging
import argparse

import pandas as pd
import pysam

from datetime import datetime

from io import StringIO
import multiprocessing as mp

LOG = logging.getLogger()
logging.basicConfig(stream=sys.stderr, level="INFO", format='%(asctime)s - %(levelname)s - %(message)s')

# These are gathered from ASD cohort
TRUE_MALEY_RATIO = 0.009847771
TRUE_FEMALEY_RATIO = 0.002462753


# datef = datetime.today().strftime('%Y%m%d')

def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    parser.add_argument("filepath", nargs=1, type=str, help="BAM file")
    parser.add_argument("--prefix", type=str, help="Prefix name will be used for output files", required=True)
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--maleY", type=int, default=TRUE_MALEY_RATIO,
                        help="The Y ratio estimated from a truely female individual")
    parser.add_argument("--femaleY", type=int, default=TRUE_FEMALEY_RATIO,
                        help="The Y ratio estimated from a truely male individual")

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    filepath = args.filepath[0]
    prefix = args.prefix
    threads = args.threads
    maleY = args.maleY
    femaleY = args.femaleY

    autosomes_set = ['chr{}'.format(x) for x in list(range(1, 23))]
    genome_dict = {entry: "autosome" for entry in autosomes_set}
    genome_dict.update({"chrY": "chrY"})

    LOG.info(f"Getting coverages in {filepath}")
    df = gather_dataframes(filepath=filepath, num_processes=threads, target_chroms=["chrY"] + autosomes_set)
    df["sample"] = filepath

    df["genome_type"] = df.apply(lambda row: genome_dict.get(row["#rname"], ""), axis=1)

    LOG.info(f"Estimating sex in {filepath}")
    stats = df.groupby(["genome_type"], group_keys=False)["numreads"].sum().reset_index().set_index("genome_type").T
    stats["ratioY"] = stats["chrY"] / stats["autosome"]
    stats["Yestimate"] = (stats["ratioY"] - femaleY) / (maleY - femaleY)
    stats["assigned_sex"] = stats["Yestimate"].apply(guesstimate_sex)
    stats["sample"] = filepath

    stats.reset_index(drop=True, inplace=True)
    df.columns.name = ""

    LOG.info(f"Done, bye")
    # os.makedirs(datef, exist_ok=True)
    # prefix = os.path.join(datef, prefix)
    stats.to_csv(f"{prefix}-Ystats.tsv", header=True, index=False, sep="\t")
    df.to_csv(f"{prefix}-stats.tsv", header=True, index=False, sep="\t")


def guesstimate_sex(x):
    if x >= 0.6:
        return "male"
    elif 0.6 > x > 0.5:
        return "maybe-male"
    else:
        return "female"


def gather_dataframes(filepath, num_processes, target_chroms):
    arglist = []
    for entry in target_chroms:
        arglist.append(make_adl_args(region=entry, min_read_len=10000))

    filepath_list = [filepath] * len(arglist)

    with mp.Pool(processes=num_processes) as pool:
        results = pool.starmap(get_coverage, zip(filepath_list, arglist))

    results = [result for result in results]
    final_result = pd.concat(results, ignore_index=True)

    return final_result


def make_adl_args(**kwargs) -> list:
    adl_args = []
    for k, v in kwargs.items():
        if "_" in k:
            k = k.replace("_", "-")
        adl_args.append(f"--{k}")
        adl_args.append(str(v))
    return adl_args


def get_coverage(filepath, adl_args: list) -> pd.DataFrame:
    try:
        out = pysam.coverage(filepath, *adl_args)
        df = pd.read_table(StringIO(out), header=0)
    except pysam.utils.SamtoolsError:
        LOG.warning(f"{adl_args} is invalid.")
        df = pd.DataFrame()
    return df


if __name__ == "__main__":
    sys.exit(main())
