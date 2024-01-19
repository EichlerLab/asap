# Housekeeping
Here you can find tasks related to SFARI deposition, quick_stats generation, and routinely updating google sheets.

##### Table of Contents
* [SFARI data deposition](#sfari)
* [Auto-populate google sheets](#google-sheets)
  * [sequencing_summary](#populate-sequencing-summary)
  * [Autism_Long_Read_Project_Yang_Mei](#populate-autism-sheets)
* [Generate quick-stats](#quick-stats)
* [FAQ](#faq)

## SFARI
Samples deposited: https://docs.google.com/spreadsheets/d/1wfmIFOv_77eGxti2wFyhMb4woBUF4M-aH0JWb2nnWtk/edit?usp=sharing
1. Please make sure the fastq.gz are cleaned of non-human reads before submission. You could do this using back-reference-qc pipeline.

Get the script [here](housekeeping_scripts/sfari_data_deposit.py) 
```shell
cd /net/eichler/vol28/projects/long_read_archive/nobackups/sharing/SFARI

# make a folder with a name corresponding to a submission batch, e.g. submission_batch_0-20230713
mkdir submission_batch_1 && cd $_

# run the script
./sfari_data_deposit --proj_dir /net/eichler/vol28/projects/long_read_archive/nobackups/clinical --sample 14694_s1:SSC12738
```

## Quick stats
The scripts are here for each [ont](housekeeping_scripts/get_ont_stats.sh) and [pacbio](housekeeping_scripts/get_pb_stats.sh)
1. Make sure you have permissions to write in the long read archive.
2. Get the above scripts and freely execute the script, e.g `./get_ont_stats.sh`
3. Set up your crontabs separately.

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
