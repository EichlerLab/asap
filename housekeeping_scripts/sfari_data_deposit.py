#!/usr/bin/env python3
"""
This is to prepare the processed SSC samples for SFARI long-read sequencing data deposition. Currently only supports PacBio + Nanopore data
Usage: sfari_data_deposit.py --sample 11918_s1:SSC --proj_dir --ont_dtype fast5
Author: Mei Wu https://github.com/projectoriented
"""
import sys

if sys.version_info < (3, 10):
    print("This script requires Python 3.10 or higher.")
    sys.exit(1)

import argparse
import os
import concurrent.futures
import pathlib
import logging
import hashlib

import pandas as pd

LOG = logging.getLogger()
logging.basicConfig(stream=sys.stdout, level="INFO", format='%(asctime)s - %(levelname)s - %(message)s')


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    parser.add_argument('--sample', metavar='14232_p1:SSC10694', required=True,
                        help='sample:ssc_name', default=None)
    parser.add_argument('--ont_dtype', required=False,
                        help='ONT raw data type', default='fast5', choices=['fast5', 'pod5'])
    parser.add_argument('--proj_dir', help='LRA directory to the samples', required=True)
    parser.add_argument(
        '--threads', help='Threads used in hard linking + generating MD5', required=False, default=8, type=int
    )
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    sample_name, ssc_name = args.sample.split(":")
    ont_raw_data_type = args.ont_dtype
    threads = args.threads

    # ------------------------------- Process ONT ------------------------------- #
    ont_obj = ONTTree(sample=sample_name, prefix=args.proj_dir, ssc_name=ssc_name, raw_data_type=ont_raw_data_type)
    ont_df = ont_obj.desired_tree_df()

    distribute_hard_link_jobs(src_list=ont_df["src_file"], dest_list=ont_df["dest_file"], threads=threads)
    ont_obj.write_out_raw_data_md5()

    ont_obj.write_bc_md5s()

    # ------------------------ Process HiFi/CCS/Revio ------------------------------- #
    hifi_obj = HiFiTree(sample=sample_name, prefix=args.proj_dir, ssc_name=ssc_name)
    hifi_df = hifi_obj.desired_tree_df()

    # William says to keep only trio hifiasm if it exists-
    trio_exists = any("trio" in x for x in hifi_df.dest_dir.unique())
    if trio_exists:
        hifi_df = hifi_df[~hifi_df.dest_dir.str.contains('hifiasm/\d+.\d+.\d+', regex=True)]

    distribute_hard_link_jobs(src_list=hifi_df["src_file"], dest_list=hifi_df["dest_file"], threads=threads)

    search_dict = hifi_obj.find_md5(df=hifi_df)
    MD5Generator(search_dict=search_dict, num_workers=threads).write_md5sums()

    dest_md5_df = get_hifi_dest_md5_df(hifi_df=hifi_df, search_dict=search_dict)

    # fix the hifiasm md5 names..
    hifiasm_df = dest_md5_df[dest_md5_df.dest_file.str.contains("hifiasm")]
    if not os.path.exists(hifiasm_df.dest_file[0]):
        file_df = pd.read_table(hifiasm_df.src_file[0], names=["md5", "file_name"])
        file_df["file_name"] = file_df.file_name.str.replace(sample_name, ssc_name)

        def manipulate_file_name(x):
            if "dip.hap1" in x:
                return x.replace("hap1", "pat")
            elif "dip.hap2" in x:
                return x.replace("hap2", "mat")
            else:
                return x

        file_df["file_name"] = file_df.file_name.apply(manipulate_file_name)
        file_df.to_csv(hifiasm_df.dest_file[0], sep="\t", index=False, header=False)
        LOG.info(f"wrote: {hifiasm_df.dest_file[0]}")
    else:
        LOG.info(f"{hifiasm_df.dest_file[0]} exists, ignoring.")

    dest_md5_df = dest_md5_df[~dest_md5_df.dest_file.str.contains("hifiasm")]

    distribute_hard_link_jobs(src_list=dest_md5_df["src_file"], dest_list=dest_md5_df["dest_file"], threads=threads)


def distribute_hard_link_jobs(src_list, dest_list, threads=8):
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit tasks for parallel execution
        futures = [executor.submit(make_dir_and_hard_link, src, dest) for src, dest in
                   zip(src_list, dest_list)]

        # Wait for all tasks to complete
        concurrent.futures.wait(futures)
        for f in futures:
            if not f.result():
                print(f)


def is_empty_dir(path_obj):
    if path_obj.is_dir():
        return not any(path_obj.iterdir())
    else:
        return True


def file_is_paired(file_list):
    num = len(file_list)
    if num % 2 == 0:
        return True
    elif num % 3 == 0:
        return True
    else:
        return False


def make_dir_and_hard_link(src_file, dest_file):
    src_path = pathlib.Path(src_file)
    dest_path = pathlib.Path(dest_file)

    dest_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        dest_path.hardlink_to(src_path)
    except FileExistsError:
        LOG.debug(f"{dest_file} exists, ignoring.")

    return 1


def get_hifi_dest_md5_df(hifi_df, search_dict):
    dest_md5_df = pd.DataFrame(columns=["src_file", "dest_file"])

    dest_md5_files=[]
    for row in hifi_df.itertuples():
        if "fastq.gz" in row.dest_file:
            dest_md5_files.append(
                os.path.join(row.dest_dir, "fastq.md5")
            )
        elif ".fasta" in row.dest_file:
            dest_md5_files.append(
                os.path.join(row.dest_dir, "fasta.md5")
            )
        elif ".bam" in row.dest_file:
            dest_md5_files.append(
                os.path.join(row.dest_dir, "bam.md5")
            )
        else:
            raise ValueError(f"hifi md5 not accnted for: {row.src_file}")
    dest_md5_files = set(dest_md5_files)

    for idx, x in enumerate(search_dict.keys()):

        # Match the last two words
        substring = os.sep.join(x.split(os.sep)[-2:])

        dest_md5 = [f for f in dest_md5_files if substring in f][0]
        dest_md5_df.loc[idx, :] = x, dest_md5

    return dest_md5_df


class MD5Generator:
    def __init__(self, search_dict, num_workers):
        self.search_dict = search_dict
        self.num_workers = num_workers

    @staticmethod
    def make_md5sum(file_name):
        md5 = hashlib.md5()
        with open(file_name, 'rb') as f:
            while True:
                chunk = f.read(4096)  # Read and hash the file in 4KB chunks
                if not chunk:
                    break
                md5.update(chunk)
        return md5.hexdigest(), os.path.basename(file_name)

    def parallize_md5sums(self, file_list):
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.num_workers) as executor:
            md5sums = list(executor.map(self.make_md5sum, file_list))
        return md5sums

    def fetch_md5sum(self):
        md5sums_dict = {}
        for k in self.search_dict.keys():
            # if not os.path.exists(k):
            fl = []
            for each_file in self.search_dict[k]:
                LOG.info(f"making md5: {each_file}")
                fl.append(self.make_md5sum(file_name=each_file))
                LOG.info(f"md5 fetched: {each_file}")
            md5sums_dict[k] = fl
            # else:
            #     LOG.debug(f"{k} exists, ignoring.")
        return md5sums_dict

    def write_md5sums(self):
        md5sums_dict = self.fetch_md5sum()
        for k in md5sums_dict.keys():
            with open(k, "a") as outfile:
                for md5, file_name in md5sums_dict[k]:
                    outfile.write(f"{md5}\t./{file_name}")
                    outfile.write("\n")
            LOG.info(f"file wrote: {k}")


class HiFiTree:
    def __init__(self, sample, prefix, ssc_name):
        self.sample = sample
        self.prefix = prefix
        self.ssc_name = ssc_name

    @staticmethod
    def get_bams(path_obj, which_one) -> list:
        first_layer = [x for x in path_obj.iterdir() if x.name.endswith('.bam') or x.name.endswith('.bam.pbi')]

        if which_one == "subread":
            bams = [x.as_posix() for x in first_layer]
        elif which_one == "hifi":
            bams = [x.as_posix() for x in first_layer if x.name.startswith('m')]
        else:
            raise ValueError(f"Unsupported which_one param: {which_one}")

        # if not file_is_paired(bams):
        #     raise OSError(f"HiFi {which_one} for {path_obj} may not all have corresponding indicies: len={len(bams)}.")

        return bams

    @staticmethod
    def get_fastq_gz(path_obj) -> list:
        first_layer = [x for x in path_obj.iterdir() if
                       x.name.endswith('.fastq.gz') or x.name.endswith('.fastq.gz.gzi') or x.name.endswith(
                           '.fastq.gz.fai')]

        fastq_gz = [x.as_posix() for x in first_layer if x.name.startswith('m')]

        # if not file_is_paired(fastq_gz):
        #     raise OSError(f"HiFi fastq.gz for {path_obj} may not all have corresponding indicies: len={len(fastq_gz)}.")

        return fastq_gz

    @property
    def get_subread_bams(self) -> list:
        path_obj = pathlib.Path(self.prefix, self.sample, "raw_data/PacBio_HiFi/subread")
        # Check if empty
        is_empty = is_empty_dir(path_obj=path_obj)
        if not is_empty:
            subread_bam_list = self.get_bams(path_obj=path_obj, which_one="subread")
        else:
            subread_bam_list = []
            LOG.info(f"No subreads for {self.sample}|{self.ssc_name}")
        return subread_bam_list

    def get_hifi_raw_data(self, file_type="bam") -> list:
        """ hifi here is a synonym for ccs,revio,hifi """
        path_obj = pathlib.Path(self.prefix, self.sample, "raw_data/PacBio_HiFi")

        # Check if empty
        is_empty = is_empty_dir(path_obj=path_obj)

        if not is_empty:
            if file_type == "bam":
                file_list = self.get_bams(path_obj=path_obj, which_one="hifi")
            elif file_type == "fastq":
                file_list = self.get_fastq_gz(path_obj=path_obj)
            else:
                raise ValueError(f"Invalid argument for file_type: {file_type}")
        else:
            file_list = []
            LOG.info(f"No hifi {file_type} for {self.sample}|{self.ssc_name}")
        return file_list

    @property
    def get_assemblies(self) -> list:
        path_obj = pathlib.Path(self.prefix, self.sample, "assemblies/hifiasm")
        # Check if empty
        is_empty = is_empty_dir(path_obj=path_obj)
        if not is_empty:
            assembly_file_list = [x.as_posix() for x in list(path_obj.rglob("*.fasta*")) if not x.is_symlink()]
            # Make sure this fix_sex_chr is not present
            assembly_file_list = [x for x in assembly_file_list if "fix_sex_chr" not in x]
        else:
            assembly_file_list = []
            LOG.info(f"No asm for {self.sample}|{self.ssc_name}")
        return assembly_file_list

    def desired_tree_df(self) -> pd.DataFrame:
        df = pd.DataFrame(
            data={
                "src_file": self.get_hifi_raw_data(file_type="fastq") + self.get_hifi_raw_data(
                    file_type="bam") + self.get_subread_bams + self.get_assemblies
            }
        )

        df["src_dir"] = df["src_file"].apply(os.path.dirname)
        for entry in df.itertuples():
            src_fn = entry.src_file

            search_word = "PacBio_HiFi"
            if "hifiasm" in src_fn:
                search_word = "hifiasm"

            target_suffix = src_fn[src_fn.index(search_word):]

            if "fastq.gz" in src_fn:
                df.loc[entry.Index, "dest_file"] = f"./fastq_assembly/{self.ssc_name}/{target_suffix}"
            elif "bam" in src_fn:
                df.loc[entry.Index, "dest_file"] = f"./raw_data/{self.ssc_name}/{target_suffix}"
            elif "hifiasm" in src_fn:
                placeholder_suffix = f"./fastq_assembly/{self.ssc_name}/genome_assembly/{target_suffix}"
                parental_origin = {"dip.hap1": "dip.pat", "dip.hap2": "dip.mat"}
                for k, v in parental_origin.items():
                    placeholder_suffix = placeholder_suffix.replace(k, v)

                placeholder_suffix = placeholder_suffix.replace(self.sample, self.ssc_name)
                df.loc[entry.Index, "dest_file"] = placeholder_suffix
            else:
                raise ValueError(f"{src_fn}")

        df["dest_dir"] = df["dest_file"].apply(os.path.dirname)
        return df

    @staticmethod
    def find_md5(df):
        md5_dict = df.groupby(["src_dir"])["src_file"].apply(list).to_dict()

        new_md5_dict = {}

        for k, v in md5_dict.items():
            if "subread" in k:
                new_key = os.path.join(k, "bam.md5")
                new_md5_dict[new_key] = md5_dict[k]

            elif "hifiasm" in k:
                new_key = os.path.join(k, "fasta.md5")
                new_md5_dict[new_key] = md5_dict[k]
            else:
                new_key = os.path.join(k, "bam.md5")
                new_md5_dict[new_key] = [x for x in md5_dict[k] if "bam" in x]
                new_key = os.path.join(k, "fastq.md5")
                new_md5_dict[new_key] = [x for x in md5_dict[k] if "fastq" in x]

        del md5_dict
        return new_md5_dict


class ONTTree:
    def __init__(self, sample, prefix, ssc_name, raw_data_type="fast5"):
        self.sample = sample
        self.prefix = prefix
        self.ssc_name = ssc_name
        self.raw_data_md5_df = pd.DataFrame()
        self.library = "STD"
        self.basecall_md5_files = []
        self.raw_data_type = raw_data_type
        self.dest_dir_list = []

    def get_data(self):
        path_obj = pathlib.Path(self.prefix, self.sample, "raw_data/nanopore", self.library)

        # Check if empty
        is_empty = is_empty_dir(path_obj=path_obj)
        if not is_empty:
            copy_record_list = [
                x.as_posix() for x in list(path_obj.rglob("*.tab.gz")) if
                "copy_record" in x.as_posix() or self.raw_data_type in x.as_posix()
            ]
        else:
            raise OSError(f"No nanopore {self.library}/{self.raw_data_type} for {self.sample}|{self.ssc_name}")

        df = pd.concat(
            [pd.read_table(x, usecols=["RUN_ID", "DEST_PATH", "SIZE", "MD5"]) for x in copy_record_list]
        ).reset_index(drop=True)

        # Only take files with raw_data_type extension.
        df = df[df.DEST_PATH.str.contains(f".{self.raw_data_type}")]

        if df.empty:
            LOG.warning(f"No files found with this suffix: .{self.raw_data_type} in {copy_record_list}. Trying ")
            sys.exit(1)
        [os.path.join(x.parent, x.name) for x in path_obj.glob(fr"*/*/{}/*.{}")]

        for entry in df.itertuples():
            path_split = entry.DEST_PATH.split(os.sep)
            target_idx = path_split.index("raw_data")

            df.loc[entry.Index, "BASECALLED_DATA_PREFIX"] = os.path.join(
                self.prefix, self.sample, *path_split[target_idx:4], "fastq", entry.RUN_ID
            )

            # Take out the first directory in case it doesn't match our sample name
            df.loc[entry.Index, "DEST_PATH"] = os.path.join(
                self.prefix, self.sample, *path_split[target_idx:]
            )

        return df

    @property
    def data_df(self):
        df = self.get_data()
        df.rename(columns={"DEST_PATH": "src_file"}, inplace=True)
        return df

    def manipulate_filename(self, file_path):
        path_obj = pathlib.Path(file_path)
        search_word = self.sample
        try:
            path_obj.name.index(self.sample)
        except ValueError:
            path_obj.name.index(self.ssc_name)
            search_word = self.ssc_name

        new_base_name = path_obj.name.replace(search_word, self.ssc_name)

        path_obj = path_obj.parent / new_base_name
        return path_obj.as_posix()

    def get_basecalled_data(self, df):
        check_these_dirs = df["BASECALLED_DATA_PREFIX"].unique().tolist()

        bc_df = pd.DataFrame()
        for entry in check_these_dirs:
            path_obj = pathlib.Path(entry)

            # Check if empty
            is_empty = is_empty_dir(path_obj=path_obj)
            if not is_empty:
                file_list = [x.as_posix() for x in list(path_obj.rglob("*")) if x.is_file() and not x.is_symlink()]
                bc_df = pd.concat([bc_df, pd.DataFrame(data={"src_file": file_list})])

        bc_df.reset_index(inplace=True, drop=True)

        md5_files = bc_df[bc_df.src_file.str.contains("\.md5")].src_file.tolist()
        self.basecall_md5_files = md5_files

        bc_df = bc_df[~bc_df.src_file.str.contains(".md5|tar.gz|.txt|.tsv.gz", regex=True)]

        for entry in bc_df.itertuples():

            src_fn = entry.src_file

            if "bam" in src_fn:
                target_suffix = src_fn[src_fn.index("fastq"):]
                target_suffix = self.manipulate_filename(target_suffix)
                bc_df.loc[entry.Index, "dest_file"] = os.path.join(".", "raw_data", self.ssc_name, "ont",
                                                                   self.library, target_suffix.replace("fastq", "bam"))
            else:
                target_suffix = src_fn[src_fn.index("ont"):]
                target_suffix = self.manipulate_filename(target_suffix)
                bc_df.loc[entry.Index, "dest_file"] = os.path.join(".", "fastq_assembly", self.ssc_name, target_suffix)

        return bc_df

    def write_bc_md5s(self):
        def get_outfile_name(df, fp, search_word):
            run_id = os.path.basename(fp).split(".")[0]
            target_dir_list = [x for x in self.dest_dir_list if run_id in x]

            if search_word == "raw_data":
                df.file_name.replace({"fastq": "bam"}, regex=True, inplace=True)
                target_name = "bam"
            else:
                target_name = "fastq"
            outfile_name = os.path.join([x for x in target_dir_list if search_word in x].pop(), f"{target_name}.md5")

            if not os.path.exists(outfile_name):
                if target_name == "fastq":
                    target_name = "fastq.gz"
                df[df.file_name.str.contains(target_name)].to_csv(outfile_name, sep="\t", header=False, index=False)
                LOG.info(f"wrote: {outfile_name}")
            else:
                LOG.debug(f"{outfile_name} exists, ignoring.")

        for file_path in self.basecall_md5_files:
            df = pd.read_csv(file_path, names=["md5", "file_name"], sep="\s+")
            if not df.file_name.str.contains(self.ssc_name).any():
                df.file_name.replace({self.sample: self.ssc_name}, regex=True, inplace=True)
            else:
                LOG.debug(f"sample name is already in its ssc form {self.sample}|ssc:{self.ssc_name}")
                # raise ValueError(f"cannot determine pattern for renaming md5 in {self.sample}: {df}")
            get_outfile_name(df=df, fp=file_path, search_word="raw_data")
            get_outfile_name(df=df, fp=file_path, search_word="fastq")

    def desired_tree_df(self):
        df = self.data_df
        bc_df = self.get_basecalled_data(df=df)
        self.dest_dir_list = bc_df["dest_file"].apply(os.path.dirname).unique().tolist()

        for entry in df.itertuples():
            src_fn = entry.src_file

            if src_fn.endswith(f".{self.raw_data_type}"):
                target_suffix = src_fn[src_fn.index("ont"):]
                df.loc[entry.Index, "dest_file"] = os.path.join(".", "raw_data", self.ssc_name, target_suffix)

        df["dest_outfile_name"] = df["dest_file"].apply(
            lambda x: os.path.join(os.path.dirname(x), f"{self.raw_data_type}.md5"))

        non_md5_related_df = pd.concat([df[["src_file", "dest_file"]], bc_df])

        df["src_basename"] = df["src_file"].apply(os.path.basename)

        def process_md5_files(row):
            return list(zip(list(row.MD5), list(row.src_basename)))

        self.raw_data_md5_df = df.groupby(
            ["dest_outfile_name"], group_keys=True
        )[["src_basename", "MD5"]] \
            .apply(process_md5_files) \
            .reset_index() \
            .rename(columns={0: "list_of_tuples"})

        return non_md5_related_df

    def write_out_raw_data_md5(self):
        for entry in self.raw_data_md5_df.itertuples():
            outfn = entry.dest_outfile_name
            if not os.path.isfile(outfn):
                with open(outfn, 'a') as outfile:
                    for t in entry.list_of_tuples:
                        md5sum = t[0]
                        file_name = f"./{t[1]}"
                        outfile.write(f"{md5sum}\t{file_name}")
                        outfile.write("\n")
                LOG.info(f"wrote {outfn} md5sums for {self.sample}|{self.ssc_name}")
            else:
                LOG.debug(f"{outfn} exists for {self.sample}|{self.ssc_name}, ignoring.")


if __name__ == "__main__":
    sys.exit(main())
