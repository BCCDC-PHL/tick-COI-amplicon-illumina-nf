#!/usr/bin/env python3
import csv
import argparse

def extract_sample_IDs(summary_file):
    extractSamples = []
    with open(summary_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['check_metrics'] == 'REVIEW' or row['check_metrics'] == 'QC_FAIL/REVIEW':
                extractSamples.append(row['sample_name'])
            else:
                pass # skip this row
        return extractSamples

def main(args):
    extractSamples = extract_sample_IDs(args.summary_file)
    with open(args.summary_outfile, 'w') as f:
        for sample in extractSamples:
            f.write(sample + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary_outfile', required=True)
    parser.add_argument('--summary_file', required=True)
    args = parser.parse_args()
    main(args)
