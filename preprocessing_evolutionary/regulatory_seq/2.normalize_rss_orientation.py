#!/usr/bin/env python3

"""
FASTA Reformatter and Classifier of RSS
This script processes FASTA files to standardize sequence orientation,
simplify sequence identifiers, and optionally classify spacer sequences
based on their length.
It supports conditional reverse-complementing of sequences depending on
filename patterns and the presence of '_REV_' tags in sequence IDs.
Spacer sequences are grouped into 12 bp, 23 bp, or other categories
and saved into separate subdirectories.

Input:
- One or more FASTA files
- Optional regex patterns to control sequence orientation handling

Output:
- Reformatted FASTA files with simplified IDs
- Spacer sequences optionally split into length-based groups (12 / 23 / other)

Author: Weronika Oprzędek
"""

from Bio import SeqIO
import argparse
import re
import os

# Regex to detect "_REV_" in sequence IDs (case-insensitive)
REV_TAG = re.compile(r"_REV_", re.IGNORECASE)

def classify_spacer(seq):
    """
    Classify spacer by its length:
    - 12: length between 9 and 15
    - 23: length between 20 and 25
    - other: all other lengths
    """
    n = len(seq)
    if 9 <= n <= 15:
        return "12"
    elif 20 <= n <= 25:
        return "23"
    else:
        return "other"
    
def shorten_id(id_str: str) -> str:
    """
    Replace long feature names in the sequence ID with shorter codes.
    """
    replacements = {
        "3'D-HEPTAMER": "3DHep",
        "5'D-HEPTAMER": "5DHep",
        "3'D-NONAMER": "3DNon",
        "5'D-NONAMER": "5DNon",
        "3'D-SPACER": "3DSpa",
        "5'D-SPACER": "5DSpa",
        "3'D-RS": "3DRS",
        "5'D-RS": "5DRS",
        "J-HEPTAMER": "JHep",
        "J-NONAMER": "JNon",
        "J-SPACER": "JSpa",
        "V-HEPTAMER": "VHep",
        "V-NONAMER": "VNon",
        "V-SPACER": "VSpa",
        "V-RS": "VRS",
        "J-RS": "JRS",
    }
    for long, short in replacements.items():
        id_str = id_str.replace(long, short)
    return id_str

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", help="FASTA files")
    ap.add_argument("--flip-if-file-matches", action="append", default=[],
                    help="regex for filenames")
    ap.add_argument("--outdir", required=True, help="output directory")
    args = ap.parse_args()

    # Compile regexes for file matching
    regexes = [re.compile(r) for r in args.flip_if_file_matches]

    for infile in args.inputs:
        fname = os.path.basename(infile)
        # Check if filename matches any of the flip regexes
        file_matches = any(r.search(fname) for r in regexes)

        out_records = []
        for rec in SeqIO.parse(infile, "fasta"):
            id_has_rev = bool(REV_TAG.search(rec.id))
            # Reverse complement if needed based on file and ID
            if (file_matches and not id_has_rev) or (not file_matches and id_has_rev):
                rec.seq = rec.seq.reverse_complement()

            # Shorten sequence ID
            rec.id = shorten_id(rec.id)
            rec.description = rec.id
            out_records.append(rec)

        # If file is a spacer, classify and group by spacer type
        if "spacer" in fname.lower():
            groups = {"12": [], "23": [], "other": []}
            for rec in out_records:
                grp = classify_spacer(str(rec.seq))
                groups[grp].append(rec)

            # Write grouped records to separate directories
            for grp, recs in groups.items():
                if not recs:
                    continue
                outpath = os.path.join(args.outdir, grp)
                os.makedirs(outpath, exist_ok=True)
                outfile = os.path.join(outpath, fname)
                SeqIO.write(recs, outfile, "fasta")
                print("saved:", outfile)
        else:
            # Write all records to output directory
            outpath = args.outdir
            os.makedirs(outpath, exist_ok=True)
            outfile = os.path.join(outpath, fname)
            SeqIO.write(out_records, outfile, "fasta")
            print("saved:", outfile)

if __name__ == "__main__":
    main()


"""
python normalize_rss_orientation.py (Get-ChildItem original/rss/*.fasta).FullName --flip-if-file-matches '^J-' --flip-if-file-matches '^5D-' --outdir fasta/rss
"""