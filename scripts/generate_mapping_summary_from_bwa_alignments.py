#!/usr/bin/env python3

import os
import subprocess
import argparse
import pandas as pd
from functools import reduce

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dirbams', type=str, default='temp/bwa/', help="Directory with BAM files to calculate statistics from")
parser.add_argument('-r', '--report', type=str, default='results/mapping_summary.csv', help='Name and path of the report file in csv format')
args = parser.parse_args()

mapping_reports = [args.dirbams + f for f in os.listdir(args.dirbams) if f.endswith('_single_mapping_summary.txt')]

mapping_report_dfs = []

for mapping_report in mapping_reports:
    sample_name = os.path.basename(mapping_report).strip('_single_mapping_summary.txt')
    temp_df = pd.read_csv(mapping_report, sep=":", header=0, names=["attribute",sample_name])
    temp_df = temp_df.replace(to_replace='\t', value='', regex=True)
    mapping_report_dfs.append(temp_df)


combined_mapping_reports = reduce(lambda left,right: pd.merge(left, right, on='attribute', how='inner'), mapping_report_dfs)
combined_mapping_reports.to_csv(path_or_buf=args.report)
