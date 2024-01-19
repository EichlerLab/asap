#!/usr/bin/env python3
"""
Usage: ./prepare.py filepath --sheet_name
Author: Mei Wu, https://github.com/projectoriented
"""
import os.path
import sys
import logging
import argparse

import pandas as pd

from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError

LOG = logging.getLogger()
logging.basicConfig(stream=sys.stdout, level="INFO", format='%(asctime)s - %(levelname)s - %(message)s')


def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    parser.add_argument('filepath', nargs=1, help='A tsv file with the following columns')
    parser.add_argument('--sheet_name', required=True, type=str, help='Sheet name in the Google Sheet')

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    filepath = args.filepath[0]
    sheet_name = args.sheet_name

    df = pd.read_table(filepath, header=0, dtype=str, keep_default_na=False)
    df.set_index(df.columns.tolist(), inplace=True)

    gs_instance = GoogleSheet()
    gs_instance.search_file(target_file="sequencing_summary")

    # Get the header of the google sub-sheet
    header_column_order = gs_instance.get_header_row(subsheet_name=sheet_name)

    # Check if there is data
    data_on_gs = gs_instance.get_values(subsheet_name=sheet_name, range_name="A2:Z")
    if data_on_gs:
        preserving_df = pd.DataFrame(data_on_gs, columns=header_column_order, dtype=str)
        target_idx = preserving_df.columns[preserving_df.columns != "SAMPLE_SUBSET"].tolist()

        preserving_df.set_index(target_idx, inplace=True)

        union_df = preserving_df.merge(df, on=preserving_df.index.names, how="outer", sort=False).fillna("")
        union_df.reset_index(inplace=True, drop=False)
        union_df = union_df[header_column_order]

        # In case there are duplicates in the union
        union_df.drop_duplicates(subset="FILE_PATH", inplace=True, keep="last")

        values = union_df.values.tolist()
        gs_instance.clear_values(range_name=f'{sheet_name}!A2:Z')
        gs_instance.populate_values(values=values, range_name=f'{sheet_name}!A2')
    else:
        df.reset_index(inplace=True, drop=False)

        df["SAMPLE_SUBSET"] = ""
        df = df[header_column_order]
        values = df.values.tolist()

        gs_instance.clear_values(range_name=f'{sheet_name}!A2:Z')
        gs_instance.populate_values(values=values, range_name=f'{sheet_name}!A2')


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


if __name__ == "__main__":
    sys.exit(main())
