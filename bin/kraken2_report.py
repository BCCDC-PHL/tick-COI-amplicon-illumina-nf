#!/usr/bin/env python3
import csv
import argparse
import os
import pandas as pd

def get_first_data_line(report,output_file,samplename,output_prefix):
    if os.path.getsize(report) > 0:
        df = pd.read_csv(report, header=None, sep = '\t') # complete kraken2 report
        df.columns = ['pct_reads_after_human_bacterial_removal', 'input_reads_after_human_bacterial_removal', 'input_reads_after_human_bacterial_removal_taxon', 'rank_code', 'ncbi_tax_id', 'scientific_name']

        df_data_line = df.iloc[[0]] # filter only unclassified info from the kraken2 report
        df_data_line.insert(0, 'run_id', output_prefix)
        df_data_line.insert(1, 'sample_id', samplename)
        df_data_line_reads = df_data_line['input_reads_after_human_bacterial_removal'].item()
        df_data_line['check_metrics'] = 'PASS' if df_data_line_reads >= 100 else 'FAIL/REPEAT'

        final_df_data_line_csv = df_data_line.to_csv(output_file, header = True, index=False)

    else:
        df = pd.DataFrame({'run_id': [output_prefix], 'sample_id':[samplename], 'pct_reads_after_human_bacterial_removal' :['NA'], 'input_reads_after_human_bacterial_removal' :['NA'], 'input_reads_after_human_bacterial_removal_taxon' :['NA'], 'rank_code' :['NA'], 'ncbi_tax_id' :['NA'], 'scientific_name':['NA']})
        df['check_metrics'] = 'FAIL/REPEAT'
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
