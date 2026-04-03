"""
Remove pseudogenes and split FASTA
This script processes a FASTA file to identify valid open reading frames (ORFs),
filter sequences to coding only, and organize them into gene-specific output files.

Input:
- FASTA file containing nucleotide sequences

Output:
- Multiple FASTA files (one per gene prefix) with cleaned sequences
- Output directory is created if it does not exist

Author: Weronika Oprzędek
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def find_best_orf(seq):
    """
    Checks 3 reading frames and selects the one without internal stop codons.
    Trims terminal stop codon if present.
    Returns the best ORF sequence or None.
    """
    seq = seq.upper().replace(" ", "").replace("\n", "")
    
    for frame in range(3):
        trimmed_seq = seq[frame:]
        trimmed_seq = trimmed_seq[:len(trimmed_seq)//3 * 3]  # trim to full codons
        if len(trimmed_seq) < 3:
            continue
        
        protein = str(Seq(trimmed_seq).translate())
        
        # check if there are no stops except possibly the last codon
        if "*" not in protein[:-1]:
            if protein.endswith("*"):
                # trim the stop codon
                trimmed_seq = trimmed_seq[:-3]
                protein = protein[:-1]
                print(f"Trimmed terminal stop codon in frame {frame+1} (length {len(seq)})")
            
            if (frame+1) != 1:
                print(f"Selected frame {frame+1} for sequence (length {len(seq)})")
            
            return trimmed_seq
    
    return None  # no frame works


def process_fasta_headers_all_frames(input_file, output_dir):
    """
    Processes a FASTA file:
    - removes pseudogenes
    - checks all 3 frames and selects ORF without stops
    - adds '_Cba' to ID
    - saves to files by prefix
    """
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_dict = {prefix: [] for prefix in ["TRBV", "TRAV", "IGHV", "IGKV", "IGLV", "TRBJ", "TRAJ", "IGHJ", "IGKJ", "IGLJ", "TRBD", "IGHD"]}
    total_input = 0
    total_kept = 0

    for record in SeqIO.parse(input_file, "fasta"):
        total_input += 1
        
        # skip pseudogenes
        if 'pseudo' in record.id.lower():
            #print(f"Skipping pseudogene: {record.id}")
            continue
        
        # find the best reading frame
        best_orf = find_best_orf(str(record.seq))
        if best_orf is None:
            print(f"Skipping {record.id}: no valid ORF in any frame")
            continue
        
        # add "_Cba"
        record.id = record.id + "_Cba"
        record.description = record.id
        
        # new record
        new_record = SeqRecord(Seq(best_orf), id=record.id, description=record.description)
        
        # assign to prefix
        prefix_found = False
        for prefix in output_dict.keys():
            if record.id.startswith(prefix):
                output_dict[prefix].append(new_record)
                prefix_found = True
                break
        
        if not prefix_found:
            print(f"WARNING: {record.id} does not match any known prefix")
        
        total_kept += 1

    # save to files
    for prefix, records in output_dict.items():
        if records:
            out_file = os.path.join(output_dir, f"{prefix}_cleaned.fasta")
            SeqIO.write(records, out_file, "fasta")
            print(f"Saved {len(records)} sequences to {out_file}")

    print(f"\nFinal statistics:")
    print(f"Processed: {total_input}")
    print(f"Kept: {total_kept}")

input_file = "camel_V.fasta"
output_dir = "processed_outputs"
process_fasta_headers_all_frames(input_file, output_dir)

input_file = "camel_D.fasta"
output_dir = "processed_outputs"
process_fasta_headers_all_frames(input_file, output_dir)

input_file = "camel_J.fasta"
output_dir = "processed_outputs"
process_fasta_headers_all_frames(input_file, output_dir)