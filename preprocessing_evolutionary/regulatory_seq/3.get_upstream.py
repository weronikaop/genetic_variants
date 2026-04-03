#!/usr/bin/env python3

"""
FASTA Normalizer and Segment Splitter of upstream immunoreceptor sequences

This script processes a FASTA file to standardize sequence orientation,
clean and normalize sequence identifiers, optionally subsample records,
and split sequences into V, D, J, and C gene segments.
It supports automatic reverse-complementing based on '_REV_' tags,
length-based filtering and trimming, and assignment of sequences to
segments based on gene name patterns. Sequences that cannot be classified
are saved separately.

Input:
- FASTA file with immunoreceptor upstream sequences

Output:
- Separate FASTA files for V, D, J, and C upstream sequences
- Additional file with unassigned sequences

Author: Weronika Oprzędek
"""

from Bio import SeqIO
import argparse
import os
import re
from collections import defaultdict
import random

# Regular expression to detect "_REV_" in IDs or descriptions (case-insensitive)
REV_TAG = re.compile(r"_REV_", re.IGNORECASE)

# Mapping of gene tokens to their segment type
TOKEN_TO_SEG = {
    "TRBV": "V", "TRAV": "V", "IGHV": "V", "IGKV": "V", "IGLV": "V",
    "TRBJ": "J", "TRAJ": "J", "IGHJ": "J", "IGKJ": "J", "IGLJ": "J",
    "TRBD": "D", "IGHD": "D",
    "TRBC": "C", "TRAC": "C", "IGH": "C", "IGKC": "C", "IGLC": "C",
}

# Create lower-cased list for fast substring checks
TOKENS_LOWER = [(tok.lower(), seg) for tok, seg in TOKEN_TO_SEG.items()]

def normalize_id(id_str: str, strip_rev=False) -> str:
    """Normalize sequence ID: optionally remove REV, strip punctuation, replace spaces with underscores."""
    if strip_rev:
        id_str = REV_TAG.sub("", id_str)
    # remove parentheses, quotes and some punctuation; replace spaces with underscore
    id_str = re.sub(r"[()'\".,;:]", "", id_str)
    id_str = re.sub(r"\s+", "_", id_str)
    id_str = re.sub(r"_+", "_", id_str)
    id_str = id_str.strip("_ ")
    return id_str

def deduplicate_ids(records):
    """Ensure all sequence IDs are unique by appending a counter if needed."""
    seen = defaultdict(int)
    for rec in records:
        seen[rec.id] += 1
        if seen[rec.id] > 1:
            new_id = f"{rec.id}_{seen[rec.id]}"
            rec.id = new_id
            rec.name = new_id
            rec.description = new_id
    return records

def process_record(rec, max_len=1000, min_len=30, strip_rev_name=False):
    """
    Process a single sequence record:
    - Flip sequence if REV is in ID or description
    - Trim from 3' end to max_len
    - Drop if too short
    - Normalize header
    """
    seq = rec.seq

    # flip if REV in id or description
    if REV_TAG.search(rec.id) or REV_TAG.search(rec.description):
        seq = seq.reverse_complement()

    # trim from 3' end (keep last max_len bases)
    if max_len is not None and len(seq) > max_len:
        seq = seq[-max_len:]

    # drop too short
    if len(seq) < min_len:
        return None

    # normalize header (optionally strip REV)
    new_id = normalize_id(rec.id, strip_rev=strip_rev_name)
    rec.id = new_id
    rec.name = new_id
    rec.description = new_id
    rec.seq = seq
    return rec

def assign_segment_by_token(name: str):
    """
    Assign segment type ('V', 'D', 'J', 'C') based on token in the name.
    Returns None if no match is found.
    """
    nl = name.upper()

    # IGH special cases
    if nl.startswith("IGHV"):
        return "V"
    if nl.startswith("IGHJ"):
        return "J"
    if re.match(r"^IGHD[0-9]", nl):
        return "D"
    if nl.startswith("IGH") and not (nl.startswith("IGHV") or nl.startswith("IGHJ") or re.match(r"^IGHD[0-9]", nl)):
        return "C"

    for tok, seg in TOKENS_LOWER:
        if tok in nl.lower():
            return seg

    return None

def split_and_write(records, outdir, min_count_for_info=False):
    """
    Split records by segment type and write to separate FASTA files.
    Also writes unassigned records to UNASSIGNED.fasta.
    Prints a summary of the split.
    """
    os.makedirs(outdir, exist_ok=True)
    writers = {seg: [] for seg in ("V","D","J","C")}
    unassigned = []

    for rec in records:
        seg = assign_segment_by_token(rec.id)
        if seg:
            writers[seg].append(rec)
        else:
            unassigned.append(rec)

    # write files
    counts = {}
    for seg, recs in writers.items():
        path = os.path.join(outdir, f"{seg}.fasta")
        if recs:
            SeqIO.write(recs, path, "fasta")
        counts[seg] = len(recs)

    # unassigned
    unp = os.path.join(outdir, "UNASSIGNED.fasta")
    if unassigned:
        SeqIO.write(unassigned, unp, "fasta")

    counts["UNASSIGNED"] = len(unassigned)

    # print summary
    print("=== Split summary ===")
    for k in ("V","D","J","C"):
        print(f"{k}: {counts[k]} sequences -> {os.path.join(outdir, k+'.fasta') if counts[k] else '(empty)'}")
    print(f"UNASSIGNED: {counts['UNASSIGNED']} sequences -> {unp if counts['UNASSIGNED'] else '(none)'}")

    # print few examples of unassigned (helpful for debugging)
    if counts["UNASSIGNED"] > 0:
        print("\nExamples of unassigned IDs (up to 20):")
        for r in unassigned[:20]:
            print(" ", r.id)
    return counts

def main():
    """
    Main entry point: parse arguments, process records, optionally subsample,
    deduplicate IDs, and split/write by segment.
    """
    p = argparse.ArgumentParser(description="Normalize FASTA, flip REV, trim from 3' end, optional random subsample, and split into V/D/J/C + UNASSIGNED.")
    p.add_argument("infile", help="input FASTA")
    p.add_argument("--outdir", default="split_by_segment", help="output directory")
    p.add_argument("--max-len", type=int, default=1000, help="max length kept (from 3' end). Default 1000")
    p.add_argument("--min-len", type=int, default=30, help="minimum length to keep. Default 30")
    p.add_argument("--strip-rev", action="store_true", help="remove REV token from header after flipping")
    p.add_argument("--sample", type=int, default=None, help="randomly select up to N sequences from the dataset")
    args = p.parse_args()

    records = []
    total = 0
    for rec in SeqIO.parse(args.infile, "fasta"):
        total += 1
        newr = process_record(rec, max_len=args.max_len, min_len=args.min_len, strip_rev_name=args.strip_rev)
        if newr:
            records.append(newr)

    print(f"Processed {total} input records -> {len(records)} kept after length filter.")

    # Random subsample if requested
    if args.sample is not None and args.sample < len(records):
        records = random.sample(records, args.sample)
        print(f"Randomly sampled {len(records)} sequences (out of {total})")

    records = deduplicate_ids(records)
    counts = split_and_write(records, args.outdir)

if __name__ == "__main__":
    main()

#python get_upstream.py original/5utr.fasta --outdir fasta --strip-rev
