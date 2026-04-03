"""
FASTA ID Cleaner for Batch Processing

This script processes all FASTA files in a specified directory,
removing parentheses from sequence identifiers and updating
descriptions to match the cleaned IDs.

It overwrites the original files with the modified records!

Input:
- Directory containing FASTA files

Output:
- Updated FASTA files with cleaned sequence IDs (in-place!)

Author: Weronika Oprzędek
"""

import os
import sys
from Bio import SeqIO

# Change working directory to the script's location
os.chdir(os.path.dirname(os.path.abspath(__file__)))
folder = sys.argv[1]  # Folder containing FASTA files, passed as a command-line argument

# Iterate over all files in the specified folder
for fname in os.listdir(folder):
    if not fname.endswith(".fasta"):
        continue  # Skip non-FASTA files
    path = os.path.join(folder, fname)

    records = []
    # Parse each FASTA file and process its records
    for rec in SeqIO.parse(path, "fasta"):
        # Remove parentheses from the record ID
        rec.id = rec.id.replace("(", "")
        rec.id = rec.id.replace(")", "")
        # Set the description to match the modified ID
        rec.description = rec.id
        records.append(rec)

    # Overwrite the original file with the modified records
    SeqIO.write(records, path, "fasta")
    print("Overwritten:", path)

# Example usage:
# python original/change_dash.py C:\Users\weronikaoprzedek\Desktop\projekt\data\noncoding/original/rss