#!/usr/bin/env python3
"""
Usage: copy-record-pod5.py pod5.md5 copy_record output_path
"""
import pandas as pd
import os
import re
import sys

input_file = sys.argv[1]
copy_record = sys.argv[2]
output_file = sys.argv[3]

actual_md5 = pd.read_table(input_file, header=None, names=["actual_md5", "pod5"])
actual_md5["pod5"]=actual_md5["pod5"].str.split(".", expand=True)[1]
other_md5=pd.read_table(copy_record, header=0)

other_md5["bn"] = other_md5.DEST_PATH.apply(os.path.basename).str.split(".", expand=True)[0]

for idx, entry in other_md5.iterrows():
    res = actual_md5.loc[actual_md5.pod5.str.contains(fr"\b{entry.bn}\b")]
    if not res.empty:
        if res.actual_md5.shape[0] > 1:
            print(res.pod5, entry.DEST_PATH)
            new_md5sum = None
        else:
            new_md5sum = res.actual_md5.tolist()[0]
        if new_md5sum:
            other_md5.loc[idx, "MD5"] = new_md5sum
            regex = re.compile(r"(?P<prefix>.*)(?P<target_dir>fast5).(?P<target_filename>.*fast5)")
            match = regex.match(entry.DEST_PATH)
            if match:
                match_dict = match.groupdict()
                match_dict["target_filename"] = match_dict["target_filename"].replace("fast5", "pod5")
                match_dict["target_dir"] = match_dict["target_dir"].replace("fast5", "pod5")
                new_path = os.path.join(*match_dict.values())
                other_md5.loc[idx, "DEST_PATH"] = new_path
            else:
                print(res.pod5, entry.DEST_PATH)

other_md5["SIZE"] = other_md5.DEST_PATH.apply(os.path.getsize)
other_md5.drop(columns=["bn"], inplace=True)

del actual_md5
parent_dir = os.path.dirname(output_file)
os.makedirs(parent_dir, exist_ok=True)
other_md5.to_csv(output_file, header=True, index=False, sep="\t")
