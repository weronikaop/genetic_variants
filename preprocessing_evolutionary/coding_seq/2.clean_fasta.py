"""
FASTA Header Processing
This script processes a FASTA file to clean and standardize headers,
filter sequences, and retain a single representative allele per gene.

Input:
- FASTA file with annotated headers (fields separated by '|')
- Required fields include gene/allele name, species, functionality, and optionally reading frame

Output:
- FASTA file with one curated sequence per gene
- Headers simplified to gene and species (e.g., TRBC1_Hsa)

Author: Weronika Oprzędek
"""

import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from collections import defaultdict

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def process_fasta_headers(input_file, output_file):
    """
    Process FASTA headers according to the specified rules:
    1. Keep only one allele for each gene
    2. Shorten species names to common abbreviations.
    3. Trim sequences to match reading frame and ensure length is divisible by 3.
    4. remove pseudogenes
    """
    gene_alleles = defaultdict(list)
    
    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            #split the header by '|'
            header_parts = record.description.split('|')
            
            if len(header_parts) < 3:
                print(f"Warning: Malformed header: {record.description}")
                continue
                
            #extract gene name (between first | and first *)
            gene_name_part = header_parts[1]
            gene_match = re.match(r'([^*]+)\*(\d+)', gene_name_part)
            
            if not gene_match:
                print(f"Warning: Could not parse gene name in header: {record.description}")
                continue
                
            gene_base = gene_match.group(1)
            allele_number = gene_match.group(2)
                
            #extract and shorten species name
            species = header_parts[2].strip()
            species = species.split()[0]
            
            #species shortening
            species_mapping = {
                'Homo': 'Hsa',
                'Mus': 'Mmu',
                'Canis': 'Cfa',
                'Sus': 'Ssc',
                'Gorilla': 'Ggo',
                'Camelus': 'Cba'
            }
            
            short_species = species_mapping.get(species, species.split()[0][0:3])

            new_header = f"{gene_base}_{short_species}"

            #remove pseudogenes
            functionality = header_parts[3].strip()
            if functionality == 'P' or functionality == '(P)':
                print(f"Skipping pseudogene: {new_header}*{allele_number}")
                continue
            
            #get reading frame - default to '1' if not specified
            reading_frame = header_parts[7].strip() if len(header_parts) > 7 else '1'
            
            sequence = str(record.seq)
            
            #trim beginning to match reading frame
            if reading_frame in ['2', '3']:
                shift = int(reading_frame) - 1
                if shift < len(sequence):
                    sequence = sequence[shift:]
                else:
                    print(f"Warning: Sequence too short for frame adjustment: {gene_base}")
                    continue
            
            #ensure sequence length is divisible by 3
            if len(sequence) % 3 != 0:
                trim_length = len(sequence) // 3 * 3
                #trim from the end if not divisible by 3
                if trim_length > 0:
                    sequence = sequence[:trim_length]
                else:
                    print(f"Skipping {gene_base}: too short after trimming ({len(sequence)} bp)")
                    continue
            
            protein = str(Seq(sequence).translate())
            if "*" in protein[:-1]:
                print(f"Skipping {new_header}*{allele_number}: internal stop codon")
                continue

            #new record
            new_record = SeqRecord(
                Seq(sequence),
                id=new_header,
                description=f"allele_{allele_number}"
            )
            
            gene_alleles[new_header].append((int(allele_number), new_record))
    
    #for each gene, select the best allele
    records_to_keep = []
    for gene_base, alleles in gene_alleles.items():
        #sort by allele number (prefer lower numbers)
        alleles.sort(key=lambda x: x[0])
        
        #select the first functional allele (lowest number)
        best_allele_number, best_record = alleles[0]
        #print(f"Selected {gene_base}*{best_allele_number:02d} for gene {gene_base}")
        
        #update the ID to include species but not allele number
        best_record.description = ""
        
        records_to_keep.append(best_record)
    
    #write to output file
    with open(output_file, "w") as output_handle:
        SeqIO.write(records_to_keep, output_handle, "fasta")
    
    print(f"Processed {len(records_to_keep)} sequences")
    print(f"Original had {sum(1 for _ in SeqIO.parse(input_file, 'fasta'))} sequences")

input_file = "TRBV.fasta"
output_file = "TRBV_headers.fasta"
process_fasta_headers(input_file, output_file)

"""#TRB
input_file = "TRBD.fasta"
output_file = "TRBD_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "TRBV.fasta"
output_file = "TRBV_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "TRBC.fasta"
output_file = "TRBC_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "TRBJ.fasta"
output_file = "TRBJ_headers.fasta"
process_fasta_headers(input_file, output_file)


# TRA
input_file = "TRAC.fasta"
output_file = "TRAC_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "TRAV.fasta"
output_file = "TRAV_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "TRAJ.fasta"
output_file = "TRAJ_headers.fasta"
process_fasta_headers(input_file, output_file)

# IGH
input_file = "IGHC.fasta"
output_file = "IGHC_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "IGHV.fasta"
output_file = "IGHV_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "IGHD.fasta"
output_file = "IGHD_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "IGHJ.fasta"
output_file = "IGHJ_headers.fasta"
process_fasta_headers(input_file, output_file)

# IGK
input_file = "IGKC.fasta"
output_file = "IGKC_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "IGKV.fasta"
output_file = "IGKV_headers.fasta"
process_fasta_headers(input_file, output_file)


# IGL
input_file = "IGLC.fasta"
output_file = "IGLC_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "IGLV.fasta"
output_file = "IGLV_headers.fasta"
process_fasta_headers(input_file, output_file)

input_file = "IGLJ.fasta"
output_file = "IGLJ_headers.fasta"
process_fasta_headers(input_file, output_file)"""