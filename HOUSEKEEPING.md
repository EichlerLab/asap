# Housekeeping
Here you can find tasks related to SFARI deposition, and routinely updating google sheets.

##### Table of Contents
* [SFARI data deposition](#sfari)
* [Auto-populate google sheets](#google-sheets)
  * [sequencing_summary](#populate-sequencing-summary)
  * [Autism_Long_Read_Project_Yang_Mei](#populate-autism-sheets)
* [Generate quick-stats](#quick-stats)
* [FAQ](#faq)

## SFARI
* Samples deposited: https://docs.google.com/spreadsheets/d/1wfmIFOv_77eGxti2wFyhMb4woBUF4M-aH0JWb2nnWtk/edit?usp=sharing
```shell
cd /net/eichler/vol28/projects/long_read_archive/nobackups/sharing/SFARI

# make a folder with a name corresponding to a submission batch, e.g. submission_batch_0-20230713
mkdir submission_batch_1 && cd $_

# run the script
./sfari_data_deposit --proj_dir /net/eichler/vol28/projects/long_read_archive/nobackups/clinical --sample 14694_s1:SSC12738
```

## Quick stats
1. Make sure you have permissions to write in the long read archive.
2. Set up your crontabs separately for each [ont](housekeeping_scripts/get_ont_stats.sh) and [pacbio](housekeeping_scripts/get_pb_stats.sh)
3. You can find the usage of the script by looking at the script.

[:arrow_double_up:](#table-of-contents)

## Google sheets
To perform any of these steps, you must first get your `credentials.json`- got here: https://developers.google.com/sheets/api/quickstart/python#enable_the_api
Any time you run a script to update the google sheets, make sure the `credentials.json` and `token.json` (this is generated the first time you do your biz) in the same working directory.

### Populate sequencing summary
[sequencing_summary](https://docs.google.com/spreadsheets/d/1zVep6eqqjfbRuvZyrpOQtYIZPCeyc2ywDwBS42BQHo8/edit#gid=0)
1. Get the entire contents of this [folder](housekeeping_scripts/gs-sequencing_summary) and the directory should look like the below.
```text
.
├── populate-gsheets_sequencing-summary.sh
├── prepare.py
└── sheets.py
```
2. To run, make sure you also have `credentials.json` in the directory before running.
```shell
export LRA=/net/eichler/vol28/projects/long_read_archive/nobackups && ./populate-gsheets_sequencing-summary.sh
```
The python scripts expects an LRA environment variable.
3. Set up your crontabs.

[:arrow_double_up:](#table-of-contents)

### Populate autism sheets
[Autism_Long_Read_Project_Yang_Mei](https://docs.google.com/spreadsheets/d/1NYBlpsY9rizPnUMT_qE1oHg8ZGnYocCaMk7CXciZqyk/edit#gid=1556958106)
As far as I know, only yangsui@uw.edu & wumei@uw.edu can write to this sheet at the moment.
1. Get the entire contents of this [folder](housekeeping_scripts/autism-sheets) and the directory should look like the below.
```text
.
├── autism-sheets.py
└── populate-autism_gsheets.sh
```
2. To run, make sure you also have `credentials.json` in the directory before running.
```shell
export LRA=/net/eichler/vol28/projects/long_read_archive/nobackups && ./populate-autism_gsheets.sh
```
The python scripts expects an LRA environment variable.
3. Set up your crontabs.

[:arrow_double_up:](#table-of-contents)

## FAQ
1. How do I check for duplicates in the google sheets?
   1. I recommend looking at the FILE_PATH column. This column should ALWAYS be unique. If there are duplicates, then check for the columns with integers or real numbers. If those are differing values between the duplicated entries then the following scenarios happened:
      1. fastq.gz files were modified (for good reason probably)
2. Are the columns with integers or numbers of string type or numeric type in Google Sheets?
   1. String type.

[:arrow_double_up:](#table-of-contents)
