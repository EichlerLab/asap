#!/usr/bin/env python3
"""
Usage: ./autism-sheets.py filepath --sheet_name
Author: Mei Wu, https://github.com/projectoriented
"""
import os.path
import sys
import logging
import argparse
import glob

import pandas as pd
import numpy as np

import pysam

from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError

LOG = logging.getLogger()
logging.basicConfig(stream=sys.stdout, level="INFO", format='%(asctime)s - %(levelname)s - %(message)s')

LRA = os.getenv("LRA")
if not LRA:
    LOG.error(f"Expecting an environment variable LRA to use for prefixing, exiting.")
    sys.exit()


def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    subparsers = parser.add_subparsers(title='commands', dest='command')

    # Prepare command
    prepare_parser = subparsers.add_parser('prepare', help='Prepare command')
    prepare_parser.add_argument('--sample_file', help='a file of sample names',
                                type=lambda fp: path_validator(fp=fp, fp_type="file"))
    prepare_parser.add_argument('--outpath', required=False, type=str, default=sys.stdout)
    prepare_parser.add_argument('--tech', required=True, type=str, choices=["PacBio_HiFi", "nanopore", "other"])
    prepare_parser.add_argument('--merqury_prefix', required=True, type=lambda fp: path_validator(fp=fp, fp_type="dir"))
    prepare_parser.add_argument('--vbi_prefix', required=True, type=lambda fp: path_validator(fp=fp, fp_type="dir"))
    prepare_parser.add_argument('--methyl_prefix', required=True, type=lambda fp: path_validator(fp=fp, fp_type="dir"))
    prepare_parser.add_argument('--illumina_fofn_prefix', required=True, type=lambda fp: path_validator(fp=fp, fp_type="dir"))
    prepare_parser.add_argument('--brqc_prefix', required=False, type=lambda fp: path_validator(fp=fp, fp_type="dir"),
                                help="back reference qc path")
    prepare_parser.set_defaults(func=prepare_cmd)

    # Sheets command
    sheets_parser = subparsers.add_parser('sheets', help='Sheets command')
    sheets_parser.add_argument('--sheet_name', help='Specify sub sheet name in Autism_Long_Read_Project_Yang_Mei')
    sheets_parser.add_argument('filepath', nargs=1, help='A tsv file with the columns in subsheet',
                               type=content_checker)
    sheets_parser.set_defaults(func=sheets_cmd)

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    args.func(args)


def sheets_cmd(args: argparse.ArgumentParser):
    filepath = args.filepath[0]
    sheet_name = args.sheet_name

    df = pd.read_table(filepath, header=0, dtype=str, keep_default_na=False)

    gs_instance = GoogleSheet()
    gs_instance.search_file(target_file="Autism_Long_Read_Project_Yang_Mei")

    # Get the header of the google sub-sheet
    header_column_order = gs_instance.get_header_row(subsheet_name=sheet_name)

    # Check if there is data
    df = df[header_column_order]
    values = df.values.tolist()

    gs_instance.clear_values(range_name=f'{sheet_name}!A2:Z')
    gs_instance.populate_values(values=values, range_name=f'{sheet_name}!A2')


def prepare_cmd(args: argparse.ArgumentParser):
    sample_file_df = pd.read_table(args.sample_file, header=0, dtype=str)
    sample_list = sample_file_df["sample"].tolist()
    tech = args.tech

    sample_file_df.rename(columns={
        "sample": "internal_id",
        "alt_id": "external_id"
    }, inplace=True)

    master_df = pd.DataFrame()
    for sn in sample_list:
        df_dict = {}
        if tech != "other":
            df_dict["quick_stats"] = combine_quick_stats(prefix=LRA, cohort="clinical", sample_name=sn, tech=tech)
            df_dict["back-reference-qc"] = get_back_reference_qc_summary_path(prefix=args.brqc_prefix, sample_name=sn,
                                                                              tech=tech)
            df_dict["vbi"] = get_vbi(prefix=args.vbi_prefix, sample_name=sn, tech=tech)
            df_dict["methyl_dss"] = get_methyl_dss_filepath(prefix=args.methyl_prefix, sample_name=sn, tech=tech)
            if tech != "nanopore":
                df_dict["merqury_stats"] = merqury_stats(prefix=args.merqury_prefix, sample_name=sn)
                df_dict["hifiasm_stats"] = hifiasm_stats(prefix=os.path.join(LRA, "clinical"), sample_name=sn)
        else:
            df_dict["other"] = get_illumina_fofn(prefix=args.illumina_fofn_prefix, sample_name=sn)

        df = pd.DataFrame()
        for v in df_dict.values():
            df = pd.concat([df, v], axis=1)

        df.fillna("", inplace=True)
        master_df = pd.concat([master_df, df])

    if tech == "other":
        sample_file_df["sample_name"] = sample_file_df["internal_id"]
        sample_file_df.set_index("sample_name", inplace=True, drop=True)
        master_df = pd.concat([master_df, sample_file_df], axis=1)

    master_df.to_csv(args.outpath, sep="\t", header=True, index=False)


def path_validator(fp, fp_type="dir"):
    if fp_type == "dir":
        if os.path.isdir(fp):
            return fp
        else:
            raise NotADirectoryError(fp)
    else:
        if os.path.exists(fp):
            return fp
        else:
            raise FileNotFoundError(fp)


def content_checker(fp):
    df = pd.read_table(fp)
    if df.empty:
        raise pd.errors.EmptyDataError()
    else:
        return fp


class GoogleCredentials:
    # If modifying these scopes, delete the file token.json.
    SCOPES = ['https://www.googleapis.com/auth/spreadsheets', "https://www.googleapis.com/auth/drive"]

    def __init__(self):
        self.creds = self.check_credentials()

    @classmethod
    def check_credentials(cls):
        creds = None
        # The file token.json stores the user's access and refresh tokens, and is
        # created automatically when the authorization flow completes for the first
        # time.
        if os.path.exists('token.json'):
            creds = Credentials.from_authorized_user_file('token.json', cls.SCOPES)
        # If there are no (valid) credentials available, let the user log in.
        if not creds or not creds.valid:
            if creds and creds.expired and creds.refresh_token:
                creds.refresh(Request())
            else:
                flow = InstalledAppFlow.from_client_secrets_file(
                    'credentials.json', cls.SCOPES)
                creds = flow.run_local_server(port=0)
            # Save the credentials for the next run
            with open('token.json', 'w') as token:
                token.write(creds.to_json())
        return creds


class GoogleSheet(GoogleCredentials):
    def __init__(self):
        super().__init__()
        self.drive_service = build('drive', 'v3', credentials=self.creds)
        self.sheet_service = build('sheets', 'v4', credentials=self.creds)
        self.spreadsheet_id = ""
        self.spreadsheet_name = ""

    def search_file(self, target_file: str) -> (str, str):
        """
        Search file in drive location
        Load pre-authorized user credentials from the environment.
        TODO(developer) - See https://developers.google.com/identity
        for guides on implementing OAuth2 for the application.
        :param target_file: Target file name
        :return: A tuple (target_file, target_file_id)
        """

        try:
            # create drive api client
            page_token = None
            file_to_return = ""

            while True:
                response = self.drive_service.files().list(
                    q="mimeType='application/vnd.google-apps.spreadsheet'",
                    spaces='drive',
                    fields='nextPageToken, ''files(id, name)',
                    pageToken=page_token
                ).execute()

                for file in response.get('files', []):
                    file_name = file.get("name")
                    file_id = file.get("id")
                    if file_name == target_file:
                        self.spreadsheet_id = file_id
                        self.spreadsheet_name = file_name
                        return file_name, file_id
                page_token = response.get('nextPageToken', None)
                if page_token is None:
                    break

        except HttpError as error:
            LOG.error(F'An error occurred: {error}')
            file_to_return = None

        return file_to_return

    def clear_values(self, range_name):
        result = self.sheet_service.spreadsheets().values().clear(
            spreadsheetId=self.spreadsheet_id,
            range=range_name
        ).execute()

        LOG.info(f"{(result.get('clearedRange'))} cells cleared.")

    def populate_values(self, values: list[list], range_name: str):

        body = {
            'values': values
        }

        result = self.sheet_service.spreadsheets().values().append(
            spreadsheetId=self.spreadsheet_id, range=range_name,
            valueInputOption='RAW', body=body).execute()
        LOG.info(f"{(result.get('updates').get('updatedCells'))} cells appended.")

    def get_values(self, subsheet_name, range_name):

        result = self.sheet_service.spreadsheets().values().get(
            spreadsheetId=self.spreadsheet_id, range=f"{subsheet_name}!{range_name}"
        ).execute()

        result = result.get('values', [])

        return result

    def get_header_row(self, subsheet_name):

        result = self.sheet_service.spreadsheets().values().get(
            spreadsheetId=self.spreadsheet_id, range=f"{subsheet_name}!A1:Z1"
        ).execute()

        result = result.get('values', [])

        if result:
            result = result.pop()

        return result


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


def combine_quick_stats(prefix, cohort, sample_name, tech="nanopore") -> pd.DataFrame:
    prefix = os.path.join(prefix, cohort, sample_name, "raw_data", tech)

    if tech == "nanopore":
        file_pattern = f'{prefix}/STD/quick_stats/n50_guppy_sup-prom*_v6.txt'
    else:
        file_pattern = f'{prefix}/quick_stats/n50_{sample_name}.tsv'

    target_columns = ["COVERAGE_X", "N50_Kbp", "FILES"]

    quick_stats = glob.glob(file_pattern)
    if len(quick_stats) == 0:
        return pd.DataFrame(
            data={'coverage_x': [np.nan], 'cell_count': [np.nan], 'internal_id': [sample_name], 'sample_source': [""]},
            index=[sample_name])

    df = pd.concat(
        [
            pd.read_table(file, header=0, engine='python',
                          dtype={"COVERAGE_X": float, "N50_Kbp": float, "CELL": str}
                          ) for file in quick_stats
        ]
    )

    # Take only the total
    df = df.query("CELL == 'total'").reset_index(drop=True)

    if df["CELL"].str.contains("total").sum() > 1:
        df = df.agg(
            {"COVERAGE_X": sum, "N50_Kbp": lambda x: sum(x) / df.shape[0], "FILES": lambda x: ",".join(x.tolist())}
        ).to_frame().T

    df = df[target_columns]

    # Make sure floats are two decimal places
    df = df.agg({
        "N50_Kbp": lambda x: f"{float(x):.2f}",
        "COVERAGE_X": lambda x: f"{float(x):.2f}",
        "FILES": lambda x: x}
    )
    df.rename(columns={"FILES": "FILE_PATH", "COVERAGE_X": "coverage_x", "N50_Kbp": "read_n50_kbp"}, inplace=True)

    for idx, row in df.iterrows():
        filepaths = row.FILE_PATH.split(",")
        df.loc[idx, "cell_count"] = len(filepaths)
        df.loc[idx, "internal_id"] = sample_name
        df.loc[idx, "sample_name"] = sample_name

        if tech == "nanopore":
            df.loc[idx, "sample_source"] = "cell-line"
        else:
            sample_src = []
            for entry in filepaths:
                bam_fp = entry.replace(".fastq.gz", ".bam")
                sam_tag = get_sam_tag(fp=bam_fp, target_tag="SM")
                if "ssc" in sam_tag.lower():
                    sample_src.append("cell-line")
                else:
                    sample_src.append("blood")

            sample_src_str = list(set(sample_src))
            sample_src_str.sort()
            df.loc[idx, "sample_source"] = "+".join(sample_src_str)

    df.drop(columns=["FILE_PATH"], inplace=True)
    df.set_index("sample_name", inplace=True, drop=True)
    return df


def hifiasm_stats(prefix, sample_name, hifiasm_ver="0.16.1"):
    phase = 'dip'
    if ('mo' in sample_name) or ('fa' in sample_name):
        phase = 'bp'

    fp = os.path.join(prefix, sample_name, "assemblies", "hifiasm", f"{sample_name}-{phase}_n50-{hifiasm_ver}.txt")

    def parse_fp(fp):
        n50_dict = {}
        with open(fp) as fp:
            for idx, line in enumerate(fp.readlines()):
                sp_l = "".join(line.split())
                if 'hap1' in sp_l:
                    hap1 = sp_l
                    n50_dict[hap1] = {}
                elif 'hap2' in sp_l:
                    hap2 = sp_l
                    n50_dict[hap2] = {}

                entries = sp_l.split(":")
                if (idx < 6) and ('hap1' not in sp_l):
                    n50_dict[hap1][entries[0]] = entries[1]
                elif (6 < idx < 13) and ('hap2' not in sp_l):
                    n50_dict[hap2][entries[0]] = entries[1]
        return n50_dict

    if not os.path.exists(fp):
        n50 = ""
    else:
        result_dict = parse_fp(fp=fp)
        d_keys = list(result_dict.keys())
        n50 = "|".join([result_dict[d_keys[0]].get('N50(Mbp)'), result_dict[d_keys[1]].get('N50(Mbp)')])

    return pd.DataFrame({f'hifiasm-{hifiasm_ver}_N50(Mbp)': n50}, index=[sample_name])


def merqury_stats(prefix, sample_name):
    # sample name workaround right now
    merqury_sn = f'{sample_name}_HIFI'

    filepaths = {}

    filepaths["qv"] = os.path.join(prefix, merqury_sn, f"{merqury_sn}.qv")
    filepaths["completeness"] = os.path.join(prefix, merqury_sn, "completeness.stats")

    for v in filepaths.values():
        try:
            path_validator(fp=v, fp_type="file")
        except FileNotFoundError:
            return pd.DataFrame({'merqury-completeness': [""], 'merqury-qv': [""]}, index=[sample_name])

    df_qv = pd.read_csv(filepaths["qv"], sep='\t', header=None, nrows=2,
                        names=['haplotype', 'assembly_kmers', 'both_kmers', 'qv', 'error_rate'])
    try:
        df_qv = df_qv.astype({'qv': int})
    except pd.errors.IntCastingNaNError:
        LOG.warning(f"{sample_name} has corrupted merqury qv output")
        return pd.DataFrame({'merqury-completeness': [""], 'merqury-qv': [""]}, index=[sample_name])

    qv = df_qv['qv'].astype(str).str.cat(sep='|')  # condense to reflect hap1|hap2
    del df_qv

    df_c = pd.read_csv(filepaths["completeness"], sep='\t', header=None,
                       names=['haplotype', 'both_kmers', 'assembly_kmers', 'read_kmers', 'completeness'])
    try:
        df_c = df_c.astype({'completeness': int})
    except pd.errors.IntCastingNaNError:
        LOG.warning(f"{sample_name} has corrupted merqury completeness output")
        return pd.DataFrame({'merqury-completeness': [""], 'merqury-qv': [""]}, index=[sample_name])

    df_c = df_c.iloc[0:2]
    c = df_c['completeness'].astype(str).str.cat(sep='|')  # condense to reflect hap1|hap2
    del df_c

    return pd.DataFrame({'merqury-completeness': c, 'merqury-qv': qv}, index=[sample_name])


def get_back_reference_qc_summary_path(prefix, sample_name, tech):
    if tech == "nanopore":
        modified_name = f"{sample_name}-ont"
    else:
        modified_name = f"{sample_name}-hifi"
    fp = os.path.join(prefix, "results", "reads_filtered", modified_name, "filtered_out", "kraken2", "summary.tsv.gz")
    if not os.path.exists(fp):
        fp = ""
    return pd.DataFrame({"back-reference-qc_summary-path": fp}, index=[sample_name])


def get_vbi(prefix, sample_name, tech):
    if tech == "nanopore":
        vbi_sn = f"{sample_name}-ont"
    elif tech == "Illumina":
        vbi_sn = f"{sample_name}-Illumina"
    else:
        vbi_sn = f"{sample_name}-hifi"

    fp = os.path.join(prefix, "vbi-db.tsv")
    df = pd.read_table(fp, header=0, names=["#SEQ_ID", "AVG_DP", "FREEMIX", "FREELK1", "FREELK0"]).rename(
        columns={"#SEQ_ID": "SEQ_ID", "FREEMIX": "verifybamid_FREEMIX"})
    df.drop_duplicates(subset=["SEQ_ID"], keep="last", inplace=True)

    df = df[["SEQ_ID", "verifybamid_FREEMIX"]]
    df = df.query(fr"SEQ_ID == '{vbi_sn}'").reset_index(drop=True)
    df.drop(columns=["SEQ_ID"], inplace=True)
    df["sample_name"] = sample_name
    df.set_index("sample_name", drop=True, inplace=True)

    if df.empty:
        df = pd.DataFrame({"verifybamid_FREEMIX": ""}, index=[sample_name])
    return df


def get_methyl_dss_filepath(prefix, sample_name, tech="nanopore"):
    tech_dict = {
        "nanopore": "ont",
        "PacBio_HiFi": "hifi"
    }

    family_id, member = sample_name.split("_")
    if ('mo' in sample_name) or ('fa' in sample_name):
        suffix = os.sep.join(["results", tech_dict[tech], "GRCh38*", "methylation", "phased", "non-trio", sample_name,
                              f"{sample_name}_cpg-pileup.hap2.bed"])
        fp = glob.glob(os.path.join(prefix, suffix))
        if len(fp) > 1:
            LOG.error(f"Found multiple methylation files for {sample_name}: {fp}")
            sys.exit()
    else:
        # return parent directory
        suffix = os.sep.join(
            ["results", tech_dict[tech], "analysis", "methylation", "GRCh38", "dss", family_id, "intersection",
             sample_name])
        fp = os.path.join(prefix, suffix)
        try:
            if len(os.listdir(fp)) <= 44:
                LOG.warning(f"{sample_name} not finished yet")
                fp = [""]
        except FileNotFoundError:
            fp = [""]

    if not fp:
        fp = [""]

    return pd.DataFrame({"methyl": fp}, index=[sample_name])


def get_illumina_fofn(prefix, sample_name):
    fp = os.path.join(prefix, f"{sample_name}_KG.fastq.fofn")
    if not os.path.exists(fp):
        fp = ""
    return pd.DataFrame({"illumina-fofn": fp}, index=[sample_name])


if __name__ == "__main__":
    # sys.argv = ["./autism-sheets.py", "--sheet_name", "db-hifi"]
    sys.exit(main())
