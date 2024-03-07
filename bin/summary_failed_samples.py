#!/usr/bin/env python3
import csv
import argparse
import os
import pandas as pd

def output_failed_samples(samplesheet, summary_report, kraken2_report, prefix):
    # Read the samplesheet file that lists all the samples included in the sequencing run
    samplesheet_ids_file = pd.read_csv(samplesheet)
    # extract sample IDs
    samples_to_check = samplesheet_ids_file["ID"].unique()

    # read the QC summary file
    qc_summary_report_file = pd.read_csv(summary_report)
    # find failed samples that are not listed in the QC summary file
    failed_samples = [sample for sample in samples_to_check if sample not in qc_summary_report_file["sample_name"].unique()]

    # read the kraken2 report
    kraken2_report_file = pd.read_csv(kraken2_report)
    # Initialize a new datafrom to store results
    results_df = pd.DataFrame(columns=kraken2_report_file.columns)
    # Iterate through the failed_samples and extract samples with matching rows from the kraken2 file
    for sample in failed_samples:
        matching_rows = kraken2_report_file[kraken2_report_file["sample_id"] == sample]
        if not matching_rows.empty:
            results_df = pd.concat([results_df, matching_rows])
            results_df.loc[results_df['check_metrics'] == 'PASS', 'reasons_failed'] = 'FAIL, fragmented assembly'
            results_df.loc[results_df['check_metrics'] == 'FAIL/REPEAT', 'reasons_failed'] = 'FAIL, no sequencing data'
        elif sample not in kraken2_report_file["sample_id"].unique():
            sample_no_data = pd.DataFrame({'sample_id': [sample]})
            results_df = pd.concat([results_df, sample_no_data], ignore_index=True)
            results_df['run_id'] = prefix
            results_df.loc[results_df['check_metrics'] != 'PASS', 'reasons_failed'] = 'FAIL, no sequencing data'


    # Write the results to a new output file
    output_file = f'{prefix}.summary.COI.failed.csv'
    results_df.to_csv(output_file, index=False)

    return results_df

def main():
    parser = argparse.ArgumentParser(description="Output a list of failed samples")
    parser.add_argument("--samplesheet", required=True, help="input sample list")
    parser.add_argument("--summary_report", required=True, help="input COI summary report")
    parser.add_argument("--kraken2_report", required=True, help="input kraken2 summary report")
    parser.add_argument("--prefix", required=True, help="outfile")
    args = parser.parse_args()

    results_df = output_failed_samples(args.samplesheet, args.summary_report, args.kraken2_report, args.prefix)

if __name__ == "__main__":
    main()

