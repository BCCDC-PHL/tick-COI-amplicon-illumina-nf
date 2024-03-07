#!/usr/bin/env python3
from Bio import Seq, SeqIO
import csv
import os
import argparse
import pandas as pd

def filter_sequences(input_file, output_file, length):
    # Parse the sequences from the input file
    sequences = list(SeqIO.parse(input_file, "fasta"))

    # Filter out sequences with the specific length
    filtered_seqs = [seq for seq in sequences if len(seq) >= length]

    # Write the filtered sequences to the output file
    output_fasta_file = SeqIO.write(filtered_seqs, output_file, "fasta")

    # Count the number of sequences in the output file
    total_count = len(filtered_seqs)
    print(f"{total_count} sequences of length >= {length} were saved in {output_fasta_file}")

    return output_file

def rename_headers(record, samplename):
    description = record.description
    new_header = f"{samplename}_{description}"
    record.id = new_header
    record.description = ""

    return record

def extract_contig_sequences(assembly_file, samplename):
    # Parse the assembly file
    assembly_records = SeqIO.parse(assembly_file, "fasta")

    # Iterate through the assembly records (contigs) and save them as individual files
    records = []
    for i, record in enumerate(assembly_records):
        contig_name = record.id
        contig_sequence = record.seq
        output_filename = f"{contig_name}_filtered.fasta"

        # Write the contig sequence to a new file
        with open(output_filename, "w") as output_file2:
            SeqIO.write(record, output_file2, "fasta")
            records.append(output_filename)

    return output_file2,records

def main():
    parser = argparse.ArgumentParser(description="Find the region with the highest coverage of snps.")
    parser.add_argument("--input_file", required=True, help="input multi-FASTA file")
    parser.add_argument("--output_file", required=True, help="output multi-FASTA file")
    parser.add_argument("--output_file_rename", required=True, help="rename header of FASTA file")
    parser.add_argument("--samplename", required=True, help="sample name")
    parser.add_argument("--length", type=int, help="length of sequence to be filtered out")
    args = parser.parse_args()

    # filter by length of sequence
    output_file = filter_sequences(args.input_file, args.output_file, args.length)

    modified_records = []

    for record in SeqIO.parse(output_file, "fasta"):
        modified_record = rename_headers(record, args.samplename)
        modified_records.append(modified_record)

    # save the modified records to a new fasta
    output_fasta_file = SeqIO.write(modified_records, args.output_file_rename, "fasta")


    # extract seqs present in the new fasta
    check_file_size = os.stat(args.output_file_rename).st_size
    if check_file_size != 0:
        output_file2,records = extract_contig_sequences(args.output_file_rename, args.samplename)
        # remove contig with the lowest read coverage
        cov_scores = [int(file.split('_cov_')[1].split('.')[0]) for file in records]
        max_score = cov_scores.index(max(cov_scores))

        for i, file in enumerate(records):
            if i != max_score:
                os.remove(file)


if __name__ == "__main__":
    main()

