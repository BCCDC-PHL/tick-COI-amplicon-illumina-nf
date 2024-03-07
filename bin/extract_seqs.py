#!/usr/bin/env python3
from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import os
import argparse
import pandas as pd

def remove_dashes_and_Ns(sequence):
    return sequence.replace('-','').replace('N', '')

def extract_contig_sequences(assembly_file, samplename):
    # Parse the assembly file
    assembly_records = SeqIO.parse(assembly_file, "fasta")

    # Iterate through the assembly records (contigs) and save them as individual files
    for i, record in enumerate(assembly_records):
        modified_seq = remove_dashes_and_Ns(str(record.seq))
        modified_record = SeqRecord(Seq(modified_seq), id=record.id, description=record.description)
        contig_name = modified_record.id
        contig_sequence = modified_record.seq
        output_filename = f"{contig_name}.consensus.fasta"

        # Write the contig sequence to a new file
        with open(output_filename, "w") as output_file2:
            SeqIO.write(modified_record, output_file2, "fasta")

    return output_file2

def main():
    parser = argparse.ArgumentParser(description="Find the region with the highest coverage of snps.")
    parser.add_argument("--alignment_file", required=True, help="alignment FASTA file")
    parser.add_argument("--samplename", required=True, help="sample name")
    args = parser.parse_args()


    # extract seqs present in the new alignment fasta file
    output_file2 = extract_contig_sequences(args.alignment_file, args.samplename)


if __name__ == "__main__":
    main()
