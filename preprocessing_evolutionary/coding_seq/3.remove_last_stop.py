"""
Terminal Stop Codon Cleaner
This script processes a FASTA file to remove terminal stop codons 
(TAA, TAG, TGA) from nucleotide sequences.

Input:
- FASTA file containing nucleotide sequences

Output:
- FASTA file with terminal stop codons removed

Author: Weronika Oprzędek
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_file = "fasta/IGHD_headers.fasta"
output_file = "fasta/IGHD_headers2.fasta"

stop_codons = {"TAA", "TAG", "TGA"}
records_to_keep = []

for record in SeqIO.parse(input_file, "fasta"):
    seq = str(record.seq).upper()
    
    #delete stop at the ends
    while len(seq) >= 3 and seq[-3:] in stop_codons:
        seq = seq[:-3]
    
    #new record
    new_record = SeqRecord(Seq(seq), id=record.id, description=record.description)
    records_to_keep.append(new_record)

#save
with open(output_file, "w") as out_handle:
    SeqIO.write(records_to_keep, out_handle, "fasta")

print(f"Processed {len(records_to_keep)} sequences. Stop codons at the end removed.")