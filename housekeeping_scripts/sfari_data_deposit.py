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
import re

import pandas as pd
pd.set_option('display.max_colwidth', None)

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
    parser.add_argument('--target', help="target data type to deposit", required=False, default="all", choices=["hifi","ont","asm","hifiont","hifiasm","all"])
    parser.add_argument('--ont_fastq_fofn', help="FASTQ FOFN for ONT.", required=False, type=str)
    parser.add_argument('--ont_bam_fofn', help="BAM FOFN for ONT.", required=False, type=str)
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    sample_name, ssc_name = args.sample.split(":")
    ont_raw_data_type = args.ont_dtype
    threads = args.threads

    # ------------------------------- Process ONT ------------------------------- #
    if args.target in ["hifiont","ont","all"]:
        if not args.ont_fastq_fofn or not args.ont_bam_fofn:
            parser.error("--ont_fastq_fofn  and --ont_bam_fofn are required when --target contains 'ont'")
            sys.exit(2)
        ont_obj = ONTTree(sample=sample_name, prefix=args.proj_dir, ssc_name=ssc_name, fastq_fofn=args.ont_fastq_fofn, bam_fofn=args.ont_bam_fofn)
        ont_df = ont_obj.desired_tree_df()

        distribute_hard_link_jobs(src_list=ont_df["src_file"], dest_list=ont_df["dest_file"], threads=threads)

        ont_bc_df = ont_df[~ont_df["src_file"].str.endswith(("fast5","pod5"))].copy()
        ont_bc_df["dest_dir"] = ont_bc_df["dest_file"].apply(os.path.dirname)
        ont_search_dict = ont_obj.find_ont_md5(df=ont_bc_df)
        
        MD5Generator(search_dict=ont_search_dict, num_workers=threads).write_md5sums()

        ont_bc_dest_md5_df = get_reads_dest_md5_df(reads_df=ont_bc_df, search_dict=ont_search_dict)
        for md5_index in range(len(ont_bc_dest_md5_df)):
            if not os.path.exists(ont_bc_dest_md5_df.dest_file[md5_index]):
                file_df = pd.read_table(ont_bc_dest_md5_df.src_file[md5_index], names=["md5", "file_name"]).sort_values(by="file_name")
                file_df["file_name"] = file_df.file_name.str.replace(sample_name, ssc_name)
                file_df.to_csv(ont_bc_dest_md5_df.dest_file[md5_index], sep="\t", index=False, header=False)
                LOG.info(f"wrote: {ont_bc_dest_md5_df.dest_file[md5_index]}")
        
        ont_obj.write_out_raw_data_md5()
        

    # ------------------------ Process HiFi/CCS/Revio ------------------------------- #

    if args.target in ["hifiont","hifi","all"]:

        hifi_obj = HiFiTree(sample=sample_name, prefix=args.proj_dir, ssc_name=ssc_name)
        hifi_df = hifi_obj.desired_tree_df()

        distribute_hard_link_jobs(src_list=hifi_df["src_file"], dest_list=hifi_df["dest_file"], threads=threads)

        hifi_search_dict = hifi_obj.find_md5(df=hifi_df)
                
        MD5Generator(search_dict=hifi_search_dict, num_workers=threads).write_md5sums()

        hifi_dest_md5_df = get_reads_dest_md5_df(reads_df=hifi_df, search_dict=hifi_search_dict)
        for md5_index in range(len(hifi_dest_md5_df)):
            if not os.path.exists(hifi_dest_md5_df.dest_file[md5_index]):
                file_df = pd.read_table(hifi_dest_md5_df.src_file[md5_index], names=["md5", "file_name"]).sort_values(by="file_name")
                file_df["file_name"] = file_df.file_name.str.replace(sample_name, ssc_name)
                file_df.to_csv(hifi_dest_md5_df.dest_file[md5_index], sep="\t", index=False, header=False)
                LOG.info(f"wrote: {hifi_dest_md5_df.dest_file[md5_index]}")
        
    # ------------------------ Process Assemblies (Hifiasm)--------------------------- #

    if args.target in ["hifiasm","asm","all"]:
        asm_obj = ASMTree(sample=sample_name, prefix=args.proj_dir, ssc_name=ssc_name)
        asm_df = asm_obj.desired_tree_df()

        
        
        # William says to keep only trio hifiasm if it exists-
        trio_exists = any("trio" in x for x in asm_df.dest_dir.unique())
        if trio_exists:
            asm_df = asm_df[~asm_df.dest_dir.str.contains('hifiasm/\d+.\d+.\d+', regex=True)]

        distribute_hard_link_jobs(src_list=asm_df["src_file"], dest_list=asm_df["dest_file"], threads=threads)

        asm_search_dict = asm_obj.find_md5(df=asm_df)
        MD5Generator(search_dict=asm_search_dict, num_workers=threads).write_md5sums()



        asm_dest_md5_df = get_asm_dest_md5_df(asm_df=asm_df, search_dict=asm_search_dict)

        # fix the hifiasm md5 names..
        hifiasm_df = asm_dest_md5_df[asm_dest_md5_df.dest_file.str.contains("hifiasm")]
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
            asm_dir = "/".join(hifiasm_df.dest_file[0].split("/")[:-1])
            if not os.path.exists(asm_dir):
                os.makedirs(asm_dir, exist_ok=True)
            file_df.to_csv(hifiasm_df.dest_file[0], sep="\t", index=False, header=False)
            LOG.info(f"wrote: {hifiasm_df.dest_file[0]}")
        else:
            LOG.info(f"{hifiasm_df.dest_file[0]} exists, ignoring.")

        distribute_hard_link_jobs(src_list=hifiasm_df["src_file"], dest_list=hifiasm_df["dest_file"], threads=threads)


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


def get_reads_dest_md5_df(reads_df, search_dict):
    dest_md5_df = pd.DataFrame(columns=["src_file", "dest_file"])

    dest_md5_files=[]
    for row in reads_df.itertuples():
        if "fastq.gz" in row.dest_file:
            dest_md5_files.append(
                os.path.join(row.dest_dir, "fastq.md5")
            )
        elif ".bam" in row.dest_file:
            dest_md5_files.append(
                os.path.join(row.dest_dir, "bam.md5")
            )
        else:
            raise ValueError(f"Reads md5 not accnted for: {row.src_file}")
    dest_md5_files = set(dest_md5_files)

    for idx, x in enumerate(search_dict.keys()):

        # Match the last two words
        substring = os.sep.join(x.split(os.sep)[-2:])

        dest_md5 = [f for f in dest_md5_files if substring in f][0]
        dest_md5_df.loc[idx, :] = x, dest_md5


    return dest_md5_df

def get_asm_dest_md5_df(asm_df, search_dict):
    dest_md5_df = pd.DataFrame(columns=["src_file", "dest_file"])

    dest_md5_files=[]
    for row in asm_df.itertuples():
        if ".fasta" in row.dest_file:
            dest_md5_files.append(
                os.path.join(row.dest_dir, "fasta.md5")
            )
        else:
            raise ValueError(f"asm md5 not accnted for: {row.src_file}")
    dest_md5_files = set(dest_md5_files)
#    print (dest_md5_files)

    for idx, x in enumerate(search_dict.keys()):

        # Match the last two words
        substring = os.sep.join(x.split(os.sep)[-2:])

        dest_md5 = [f for f in dest_md5_files if substring in f][0]
        dest_md5_df.loc[idx, :] = x, dest_md5


    return dest_md5_df


def get_md5sum(file_name):
    md5 = hashlib.md5()
    with open(file_name, 'rb') as f:
        while True:
            chunk = f.read(4096)  # Read and hash the file in 4KB chunks
            if not chunk:
                break
            md5.update(chunk)
    return md5.hexdigest()


class MD5Generator:
    def __init__(self, search_dict, num_workers=1):
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

    def fetch_md5sum(self):
        def get_rerun_flag(search_dict):
            rerun_flag = 0
            for fofn in search_dict.keys():
                if not os.path.isfile(fofn):
                    rerun_flag += 1
                    continue
                target_list = search_dict[fofn]
                base_names = [ file_path.split("/")[-1] for file_path in target_list ]
                fofn_df = pd.read_csv(fofn, sep="\t",header=None, names=["md5","fname"])
                fofn_df["fname"] = fofn_df["fname"].str.replace("^./","", regex=True)
                fofn_df["base_name"] = fofn_df["fname"].apply(os.path.basename)
                fofn_df["base_name"] = fofn_df["base_name"].str.replace("^./","", regex=True)

                for fname in fofn_df["fname"]:
                    if not fname in base_names:
                        rerun_flag += 1
                fofn_df_base_names = fofn_df["base_name"].tolist()
                for base_name in base_names:
                    if not base_name in fofn_df_base_names:
                        rerun_flag += 1
            return rerun_flag


        md5sums_dict = {}
        rerun_flag = get_rerun_flag(self.search_dict)
        for k in self.search_dict.keys():
            # if not os.path.exists(k):
            fl = []
            if isinstance(self.search_dict[k], list):
                search_space = self.search_dict[k]
            else:
                search_space = [self.search_dict[k]]
            if rerun_flag > 0:
                for each_file in search_space:
                    LOG.info(f"making md5: {each_file}")
                    fl.append(self.make_md5sum(file_name=each_file))
                    LOG.info(f"md5 fetched: {each_file}")
                
            else:
                LOG.info("md5 exists, skipping calculation.")
                exist_md5sum_df = pd.read_csv(k, sep="\t", header=None).drop_duplicates().sort_values(by=1)
                exist_md5sum_df[1] = exist_md5sum_df[1].str.replace("^./","",regex=True)
                fl = list(exist_md5sum_df.itertuples(index=False, name=None))
            md5sums_dict[k] = fl
                
        return md5sums_dict

    def write_md5sums(self):
        md5sums_dict = self.fetch_md5sum()
        for k in md5sums_dict.keys():
            with open(k, "w") as outfile:
                for md5, file_name in md5sums_dict[k]:
                    outfile.write(f"{md5}\t./{file_name}")
                    outfile.write("\n")
            LOG.info(f"file wrote: {k}")

class ASMTree:
    def __init__(self, sample, prefix, ssc_name):
        self.sample = sample
        self.prefix = prefix
        self.ssc_name = ssc_name

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
                "src_file": self.get_assemblies
            }
        )

        df["src_dir"] = df["src_file"].apply(os.path.dirname)
        for entry in df.itertuples():
            src_fn = entry.src_file

            search_word = "hifiasm"
            target_suffix = src_fn[src_fn.index(search_word):]

            if "hifiasm" in src_fn:
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
            if "hifiasm" in k:
                new_key = os.path.join(k, "fasta.md5")
                new_md5_dict[new_key] = md5_dict[k]

        del md5_dict
        return new_md5_dict


class HiFiTree:
    def __init__(self, sample, prefix, ssc_name):
        self.sample = sample
        self.prefix = prefix
        self.ssc_name = ssc_name

    @staticmethod
    def get_bams(path_obj, which_one) -> list:
        first_layer = [x for x in path_obj.iterdir() if (x.name.endswith('.bam') or x.name.endswith('.bam.pbi')) and (not "fail" in x.name)]

        if which_one == "hifi":
            bams = [x.as_posix() for x in first_layer if x.name.startswith('m')]
        else:
            raise ValueError(f"Unsupported which_one param: {which_one}")

        return bams

    @staticmethod
    def get_fastq_gz(path_obj) -> list:
        first_layer = [x for x in path_obj.iterdir() if
                       x.name.endswith('.fastq.gz') or x.name.endswith('.fastq.gz.gzi') or x.name.endswith('.fastq.gz.fai') and (not "fail" in x.name)]

        fastq_gz = [x.as_posix() for x in first_layer if x.name.startswith('m')]

        return fastq_gz

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

    def desired_tree_df(self) -> pd.DataFrame:

        df = pd.DataFrame(
            data={
                "src_file": self.get_hifi_raw_data(file_type="fastq") + self.get_hifi_raw_data(
                    file_type="bam")
            }
        )

        df["src_dir"] = df["src_file"].apply(os.path.dirname)
        for entry in df.itertuples():
            src_fn = entry.src_file

            search_word = "PacBio_HiFi"
            target_suffix = src_fn[src_fn.index(search_word):]

            if "fastq.gz" in src_fn:
                df.loc[entry.Index, "dest_file"] = f"./fastq_assembly/{self.ssc_name}/{target_suffix}"
            elif "bam" in src_fn:
                df.loc[entry.Index, "dest_file"] = f"./raw_data/{self.ssc_name}/{target_suffix}"
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
            else:
                new_key = os.path.join(k, "bam.md5")
                new_md5_dict[new_key] = [x for x in md5_dict[k] if "bam" in x]
                new_key = os.path.join(k, "fastq.md5")
                new_md5_dict[new_key] = [x for x in md5_dict[k] if "fastq" in x]

        del md5_dict
        return new_md5_dict


class ONTTree:
    def __init__(self, sample, prefix, ssc_name, fastq_fofn, bam_fofn):
        self.sample = sample
        self.prefix = prefix
        self.ssc_name = ssc_name
        self.raw_data_md5_df = pd.DataFrame()
        self.library = "STD"
        self.basecall_md5_files = []
        self.dest_dir_list = []
        self.fastq_fofn = fastq_fofn
        self.bam_fofn = bam_fofn


    def get_data(self):

        def sort_cr_list(copy_record_list):

            sorted_copy_record_list = [ copy_record for copy_record in copy_record_list if "/pod5/" in copy_record ]
            pod5_run_ids = [ copy_record.split("/")[-1].split(".")[0] for copy_record in sorted_copy_record_list ]
            fast5_copy_record_list = [ copy_record for copy_record in copy_record_list if "/fast5/" in copy_record ]

            for copy_record in fast5_copy_record_list:
                run_id = copy_record.split("/")[-1].split(".")[0]
                if not run_id in pod5_run_ids:
                    sorted_copy_record_list.append(copy_record)
            return sorted_copy_record_list            


        path_obj = pathlib.Path(self.prefix, self.sample, "raw_data/nanopore", self.library)


        # Check if empty
        is_empty = is_empty_dir(path_obj=path_obj)
        df = pd.DataFrame()
        if not is_empty:
            search_for_cr = list(path_obj.rglob(fr"copy_record/*/*.tab.gz"))
            if search_for_cr:
                copy_record_list = sort_cr_list([x.as_posix() for x in search_for_cr])
                df = pd.concat(
                    [pd.read_table(x, usecols=["RUN_ID", "DEST_PATH", "SIZE", "MD5"], dtype='object') for x in
                     copy_record_list]
                ).reset_index(drop=True)

                for fp in copy_record_list:
                    check_df = pd.read_table(fp, usecols=["RUN_ID", "DEST_PATH", "SIZE", "MD5"], dtype='object')
                    for x in check_df.itertuples():
                        search_path = os.path.join(self.prefix, x.DEST_PATH)
                        exists = os.path.exists(search_path)
                        if not exists:
                            raise OSError(f"{search_path} doesn't exist.")
            else:
                LOG.warning(
                    f"No copy record found. Revert to generating MD5 sums."
                )
                md5_these_filepaths = path_obj.glob(fr"*/*/*/*.{fast5,pod5}")
                for idx, fp in enumerate(md5_these_filepaths):
                    raw_data_type = fp.split(".")[-1]
                    runid = fp.absolute().parts[13]
                    parent_path_len = len(path_obj.parts) - 4
                    dest_path = os.path.join(*fp.parts[parent_path_len:])
                    size = os.path.getsize(fp)
                    md5sum_dict = MD5Generator(
                        search_dict={raw_data_type: os.path.join(fp.parent, fp.name)}).fetch_md5sum()
                    md5sum = md5sum_dict[raw_data_type][0][0]

                    # Populate the DF
                    df.loc[idx, ["RUN_ID", "DEST_PATH", "SIZE", "MD5"]] = runid, dest_path, size, md5sum

        else:
            raise OSError(f"No nanopore {self.library} for {self.sample}|{self.ssc_name}")

        # Only take files with raw_data_type extension.

        df["DEST_PATH"] = df["DEST_PATH"].fillna("")
        df['DEST_PATH'] = df['DEST_PATH'].str.replace(f"{self.prefix}/", "")
        df['DEST_PATH'] = self.prefix + "/" + df['DEST_PATH']
        df = df[df.DEST_PATH.str.contains(r"\.(?:fast5|pod5)$", regex=True)]

        for entry in df.itertuples():
            path_split = entry.DEST_PATH.split(os.sep)
            target_idx = path_split.index("raw_data")


            df.loc[entry.Index, "BASECALLED_DATA_PREFIX"] = os.path.join(
                self.prefix, self.sample, *path_split[target_idx:target_idx+3], "fastq", entry.RUN_ID
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
    
    @staticmethod
    def find_ont_md5(df):
        md5_dict = df.groupby(["src_dir"])["src_file"].apply(list).to_dict()
        new_md5_dict = {}

        for basecall_dir, basecalled_files in md5_dict.items():
            fastq_files = []
            bam_files = []
            for bc_file in basecalled_files:
                if re.search("fastq.gz", bc_file):
                    fastq_files.append(bc_file)
                elif re.search(".bam", bc_file):
                    bam_files.append(bc_file)

            fastq_key = os.path.join(basecall_dir, "fastq.md5")
            bam_key = os.path.join(basecall_dir, "bam.md5")
            new_md5_dict[fastq_key] = fastq_files
            new_md5_dict[bam_key] = bam_files

        del md5_dict
        return new_md5_dict
    
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
        

        md5_files = bc_df[bc_df.src_file.str.contains(".md5$")].src_file.tolist()

        self.basecall_md5_files = md5_files

        # bc_df = bc_df[~bc_df.src_file.str.contains(".md5\b|tar.gz\b|.txt\b|.tsv.gz\b|.html\b", regex=True)]
        selected_fastq_list = pd.read_csv(self.fastq_fofn, header=None)[0].tolist()
        selected_bam_list = pd.read_csv(self.bam_fofn, header=None)[0].tolist()
        passed_basecalled_data_list = []
        
        for fastq in selected_fastq_list:
            passed_basecalled_data_list.append(fastq)
            passed_basecalled_data_list.append(fastq+".fai")
            passed_basecalled_data_list.append(fastq+".gzi")
        for bam in selected_bam_list:
            passed_basecalled_data_list.append(bam)
            passed_basecalled_data_list.append(bam+".bai")        

        bc_df = bc_df[
            bc_df.src_file.isin(passed_basecalled_data_list)
            ]
        
        bc_df["src_dir"] = bc_df["src_file"].apply(os.path.dirname)


        for entry in bc_df.itertuples():

            src_fn = entry.src_file

            if ".bam" in src_fn: # BAM FILE
                target_suffix = src_fn[src_fn.index("fastq"):]
                target_suffix = self.manipulate_filename(target_suffix)
                bc_df.loc[entry.Index, "dest_file"] = os.path.join(".", "raw_data", self.ssc_name, "nanopore", self.library, target_suffix.replace("fastq", "bam"))
            else:
                target_suffix = src_fn[src_fn.index("nanopore"):]
                target_suffix = self.manipulate_filename(target_suffix)
                bc_df.loc[entry.Index, "dest_file"] = os.path.join(".", "fastq_assembly", self.ssc_name, target_suffix)

        return bc_df

    def desired_tree_df(self):

        def define_out_dir(dest_file):
            if "/pod5/" in dest_file:
                raw_data_type = "pod5"
            else:
                raw_data_type = "fast5"
            return os.path.join(os.path.dirname(dest_file), f"{raw_data_type}.md5")

        df = self.data_df
        bc_df = self.get_basecalled_data(df=df)
        self.dest_dir_list = bc_df["dest_file"].apply(os.path.dirname).unique().tolist()

        for entry in df.itertuples():
            src_fn = entry.src_file


            target_suffix = src_fn[src_fn.index("nanopore"):]
            df.loc[entry.Index, "dest_file"] = os.path.join(".", "raw_data", self.ssc_name, target_suffix)

        df["dest_outfile_name"] = df["dest_file"].apply(define_out_dir)
        
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
