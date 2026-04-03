"""
Exon Combiner for Immunoreceptor C regions

This script processes a FASTA file containing exon sequences and combines
exons belonging to the same gene into a single contiguous sequence.

Input:
- FASTA file where each record represents an exon
- Headers must contain allele information in the second '|' separated field
  (e.g., >...|TRBV1*01|...)

Output:
- FASTA file with one combined sequence per allele

Author: Weronika Oprzędek
"""

from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def combine_exons(input_file, output_file):
    """
    Combine exons for each allele in the order they appear in the file
    Keep the original header from the first exon
    """
    allele_exons = defaultdict(list)  # Store exons and headers for each allele
    
    # Read and group exons by allele
    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            header_parts = record.description.split('|')
            
            if len(header_parts) < 2:
                continue
                
            # Extract allele information from second field
            allele_name = header_parts[1]  # TRBC1*01, TRBV1*01, etc.
            
            # Clean sequence of newlines and spaces
            sequence = str(record.seq).replace('\n', '').replace(' ', '')
            
            # Store both sequence and the full header
            allele_exons[allele_name].append((record.description, sequence))
    
    # Combine exons for each allele
    combined_records = []
    
    for allele_name, exon_data in allele_exons.items():
        # Combine sequences in the order they appear in the file
        combined_seq = ''.join([seq for _, seq in exon_data])
        
        # Use the header from the first exon
        first_header = exon_data[0][0]
        
        # Create new record with original header
        new_record = SeqRecord(
            Seq(combined_seq),
            id=first_header,  # Keep the original header
            description=f"combined_{len(exon_data)}_exons"
        )
        
        combined_records.append(new_record)
        print(f"Combined {allele_name}: {len(exon_data)} exons → {len(combined_seq)} nt")
    
    # Write to file
    with open(output_file, "w") as output_handle:
        SeqIO.write(combined_records, output_handle, "fasta")
    
    print(f"Created {len(combined_records)} combined sequences")

combine_exons("human/trbc.fasta", "human/trbc_combined.fasta")
combine_exons("human/trac.fasta", "human/trac_combined.fasta")
combine_exons("human/ighc.fasta", "human/ighc_combined.fasta")
combine_exons("human/iglc.fasta", "human/iglc_combined.fasta")
combine_exons("human/igkc.fasta", "human/igkc_combined.fasta")

combine_exons("dog/trbc.fasta", "dog/trbc_combined.fasta")
combine_exons("dog/trac.fasta", "dog/trac_combined.fasta")
combine_exons("dog/ighc.fasta", "dog/ighc_combined.fasta")
combine_exons("dog/iglc.fasta", "dog/iglc_combined.fasta")
combine_exons("dog/igkc.fasta", "dog/igkc_combined.fasta")

combine_exons("mouse/trbc.fasta", "mouse/trbc_combined.fasta")
combine_exons("mouse/trac.fasta", "mouse/trac_combined.fasta")
combine_exons("mouse/ighc.fasta", "mouse/ighc_combined.fasta")
combine_exons("mouse/iglc.fasta", "mouse/iglc_combined.fasta")
combine_exons("mouse/igkc.fasta", "mouse/igkc_combined.fasta")

combine_exons("pig/trbc.fasta", "pig/trbc_combined.fasta")
combine_exons("pig/ighc.fasta", "pig/ighc_combined.fasta")
combine_exons("pig/iglc.fasta", "pig/iglc_combined.fasta")
combine_exons("pig/igkc.fasta", "pig/igkc_combined.fasta")

combine_exons("gorilla/trbc.fasta", "gorilla/trbc_combined.fasta")
combine_exons("gorilla/trac.fasta", "gorilla/trac_combined.fasta")
combine_exons("gorilla/ighc.fasta", "gorilla/ighc_combined.fasta")
combine_exons("gorilla/iglc.fasta", "gorilla/iglc_combined.fasta")
combine_exons("gorilla/igkc.fasta", "gorilla/igkc_combined.fasta")