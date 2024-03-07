#!/usr/bin/env python3

import csv
import argparse
import os

def get_first_data_line(blastn):
    with open(blastn, newline='') as csvfile:
        reader = csv.reader(csvfile)
        # Get the header line
        header = next(reader)
        # Get the first data line
        try:
            data_line = next(reader)
        except StopIteration:
            data_line = None
        return data_line

def main(args):
    new_header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    top_hits = []
    for blastn in args.blastn:
        if os.path.getsize(blastn) > 0:
            blastn_top_hit = get_first_data_line(blastn)
            if blastn_top_hit is not None:
                top_hits.append(blastn_top_hit)
            else:
                top_hits.append([])
        else:
            top_hits.append([])
    with open(args.outfile, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(new_header)
        writer.writerows(top_hits)
    
    header_dict = {
        'qseqid': 'Query sequence ID (generated amplicon consensus sequence)',
        'sseqid': 'Subject sequence ID (reference sequence ID in database)',
        'pident': 'Percent of identical matches',
        'length': 'Alignment length',
        'mismatch': 'Number of mismatches',
        'gapopen': 'Number of gap openings',
        'qstart': 'Start of alignment in query sequence',
        'qend': 'End of alignment in query sequence',
        'sstart': 'Start of alignment in subject sequence',
        'send': 'End of alignment in subject sequence',
        'evalue': 'Expect value or significance threshold for the match',
        'bitscore': 'Bit score or normalized score for the match'
    }
    with open(args.dictionary, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['Header label', 'Meaning'])
        writer.writeheader()
        for key, value in header_dict.items():
            writer.writerow({'Header label': key, 'Meaning': value})

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--dictionary', required=True)
    parser.add_argument('--blastn', required=True, nargs='+')
    args = parser.parse_args()
    main(args)
