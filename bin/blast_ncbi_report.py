#!/usr/bin/env python3
import csv
import argparse
import os
import pandas as pd

def get_first_data_line(report,output_file,samplename,output_prefix):
    if os.path.getsize(report) > 0:
        df = pd.read_csv(report, header=None, sep = '\t') # complete report
        df.columns = ['qseqid', 'staxids','sseqid', 'sscinames', 'scomnames', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

        df_data_line = df.iloc[[0]] # filter only top result
        df_data_line.insert(0, 'run_id', output_prefix)
        df_data_line.insert(1, 'sample_id', samplename)
        df_data_line_reads = float(df_data_line['pident'].item())
        df_data_line_qend = df_data_line['qend']
        df_data_line['amplicon_size'] = df_data_line_qend
        df_data_line['check_metrics'] = 'PASS' if df_data_line_reads >= 97 else 'QC_FAIL/REVIEW'
        final_df_data_line_csv = df_data_line.to_csv(output_file, header = True, index=False)

    else:
        df = pd.DataFrame({'run_id': [output_prefix], 'sample_id':[samplename], 'qseqid' :['NA'], 'staxids':['NA'], 'sseqid' :['NA'], 'sscinames':['NA'], 'scomnames':['NA'], 'pident' :['NA'], 'length' :['NA'], 'mismatch' :['NA'], 'gapopen':['NA'], 'qstart':['NA'], 'qend':['NA'], 'sstart':['NA'], 'send':['NA'], 'evalue':['NA'], 'bitscore':['NA']})
        df['check_metrics'] = 'FAIL'
        final_df_csv = df.to_csv(output_file, header = True, index=False)

    return output_file

def main():
    parser = argparse.ArgumentParser(description="Find the region with the highest coverage of snps.")
    parser.add_argument("--report", required=True, help="input multi-FASTA file")
    parser.add_argument("--output_file", required=True, help="output report name")
    parser.add_argument("--samplename", required=True, help="sample name")
    parser.add_argument("--prefix", required=True, help="sequencing run ID")
    args = parser.parse_args()

    output_file = get_first_data_line(args.report, args.output_file, args.samplename, args.prefix)

if __name__ == "__main__":
    main()
