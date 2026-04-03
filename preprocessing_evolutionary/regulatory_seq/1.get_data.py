"""
EMBL Feature Extractor
This script parses EMBL files to extract regulatory and structural features 
and saves them into grouped FASTA files.

Input:
- One or more EMBL files with IMGT-style annotations

Output:
- Multiple FASTA files (one per feature category) in the output directory
- Optional text file listing all feature types present in an EMBL file

Author: Weronika Oprzędek
"""

from Bio import SeqIO
import os
import warnings

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def list_feature_types(embl_file, out_txt="feature_types.txt"):
    record = SeqIO.read(embl_file, "embl")
    feature_types = sorted(set(f.type for f in record.features))

    with open(out_txt, "w") as out:
        for ftype in feature_types:
            out.write(ftype + "\n")

    print(f"[OK] Saved {len(feature_types)} feature types to {out_txt}")
    return feature_types


species_map = {
    "Homo sapiens (human)": "Hsa",
    "Mus musculus (house mouse)": "Mmu",
    "Gorilla gorilla gorilla (western lowland gorilla)": "Ggo",
    "Canis lupus familiaris (dog)": "Cfa",
    "Sus scrofa (pig)": "Ssc",
    "Camelus dromedarius": "Cdr"
}

# regulatory features and corresponding file names
regulatory_features = {
    "5'UTR": "5UTR",
    #"L-V-GENE-UNIT": "LVUNIT",
    "L-INTRON-L": "L-INTRON-L",
    "V-INTRON": "V-INTRON",
    "INTRON": "INTRON",
    "V-HEPTAMER": "V-HEPTAMER",
    "V-NONAMER": "V-NONAMER",
    "V-SPACER": "V-SPACER",
    "V-RS": "V-RSS",
    "3'D-HEPTAMER": "3D-HEPTAMER",
    "3'D-NONAMER": "3D-NONAMER",
    "3'D-SPACER": "3D-SPACER",
    "3'D-RS": "3D-RSS",
    "5'D-HEPTAMER": "5D-HEPTAMER",
    "5'D-NONAMER": "5D-NONAMER",
    "5'D-SPACER": "5D-SPACER",
    "5'D-RS": "5D-RSS",
    "J-HEPTAMER": "J-HEPTAMER",
    "J-NONAMER": "J-NONAMER",
    "J-SPACER": "J-SPACER",
    "J-RS": "J-RSS",
    "ACCEPTOR-SPLICE": "ACCEPTOR-SPLICE",
    "DONOR-SPLICE": "DONOR-SPLICE",
    "INIT-CODON": "INIT",
    "J-C-INTRON": "J-C-INTRON",
    #"POLYA_SIGNAL": "POLYA",
    #"POLYA_SITE": "POLYA",
    "INT-DONOR-SPLICE": "DONOR-SPLICE",
    "OCTAMER": "OCTAMER"
}

#list of known types
known_features = [
    '1st-CYS','2nd-CYS',"3'D-HEPTAMER","3'D-NONAMER","3'D-RS","3'D-SPACER","3'UTR","3'V-REGION",
    "5'D-HEPTAMER","5'D-NONAMER","5'D-RS","5'D-SPACER","5'J-REGION","5'UTR",
    'ACCEPTOR-SPLICE','C-CLUSTER','C-GENE','C-GENE-UNIT','CDR1-IMGT','CDR2-IMGT','CDR3-IMGT',
    'CO-PART1','CO-PART2','CONSERVED-TRP','CYTOPLASMIC-REGION','D-CLUSTER','D-GENE',
    'D-GENE-UNIT','D-J-C-CLUSTER','D-J-CLUSTER','D-REGION','DONOR-SPLICE',
    'EX0','EX1','EX2','EX3','EX4','FR1-IMGT','FR2-IMGT','FR3-IMGT',
    'IMGT-LOCUS-UNIT','INIT-CODON','INTRON','J-C-CLUSTER','J-C-INTRON','J-CLUSTER',
    'J-GENE','J-GENE-UNIT','J-HEPTAMER','J-MOTIF','J-NONAMER','J-PHE','J-REGION','J-RS','J-SPACER',
    'L-INTRON-L','L-PART1','L-PART2','L-V-GENE-UNIT','STOP-CODON','TM-PART1','TM-PART2',
    'V-CLUSTER','V-D-J-C-CLUSTER','V-EXON','V-GENE','V-HEPTAMER','V-INTRON','V-NONAMER',
    'V-REGION','V-RS','V-SPACER', 'GAP', 'CH1', 'CH2', 'CH3', 'CH4', 'CH4-CHS', 'CHS', 'M1', 'M2',
    'MISC_FEATURE','REPEAT_UNIT','DECAMER','VARIATION', 'INT-DONOR-SPLICE', 'CH3-CHS',
    'GENE','GENE-UNIT','V-D-CLUSTER','V-D-J-CLUSTER', 'V-J-C-CLUSTER',
    'POLYA_SIGNAL','POLYA_SITE', 'TRANSMEMBRANE-RE', 'EX4UTR', 'J-TRP', 'CONNECTING-REGION',
    'H1', 'H2', 'H3','H4', 'H', 'H-CH2', 'M', 'L-GENE', 'L-GENE-UNIT', 'V-LIKE', 'BC-LOOP', 'C\'C\"-LOOP',
    'FG-STRAND', 'CL', 'C-REGION', 'CDR3', 'INSERTION', 'V-J-GENE', 'L-V-J-GENE-UNIT', 'V-J-EXON', 'V-J-REGION',
    'JUNCTION', 'N-REGION', 'FR4-IMGT', 'OCTAMER'
    ]

def parse_multiple_embl(embl_files, out_dir="out_fasta"):
    os.makedirs(out_dir, exist_ok=True)
    collected = {}

    for embl_file in embl_files:
        record = SeqIO.read(embl_file, "embl")
        species = record.annotations.get("organism", "Unknown")
        species_code = species_map.get(species, species[:3])

        last_gene = None  #remember the last known gene

        for feature in record.features:
            ftype = feature.type

            #determine gene name
            gene_name = None
            for qual in ["IMGT_gene", "gene"]:
                if qual in feature.qualifiers:
                    gene_name = feature.qualifiers[qual][0]
                    last_gene = gene_name
                    break
            if not gene_name and last_gene:
                gene_name = last_gene
            if not gene_name:
                gene_name = "UnknownGene"

            if not gene_name.startswith("IGL") and not gene_name.startswith("IGK") and not gene_name.startswith("IGH") and not gene_name.startswith("TRA") and not gene_name.startswith("TRB"):
                print(f"Skipping non-TCR/BCR gene: {gene_name} in {embl_file}")
                continue

            #regulatory → save
            if ftype in regulatory_features:
                group = regulatory_features[ftype]
                seq = feature.extract(record.seq)

                #reverse complement for strand -
                if feature.location.strand == -1:
                    seq = seq.reverse_complement()
                    header = f">{gene_name}_{ftype}_REV_{species_code}"
                else:
                    header = f">{gene_name}_{ftype}_{species_code}"

                collected.setdefault(group, []).append(f"{header}\n{seq}")

            #unknown feature → warning
            elif ftype not in known_features:
                warnings.warn(f"[{species_code}] Unknown feature: {ftype} (file: {embl_file})")

        #eof

    # save to separate FASTA files
    for group, seqs in collected.items():
        fname = os.path.join(out_dir, f"{group}.fasta")
        with open(fname, "w") as out:
            out.write("\n".join(seqs))
        print(f"[OK] saved {len(seqs)} sequences to {fname}")


"""file_list = ["imgt/IGL_gorilla.txt", "imgt/IGL_human.txt", "imgt/IGL_mouse.txt","imgt/IGL_dog.txt", "imgt/IGL_pig.txt",
             "imgt/IGK_gorilla.txt", "imgt/IGK_human.txt", "imgt/IGK_mouse.txt","imgt/IGK_dog.txt", "imgt/IGK_pig.txt",
             "imgt/IGH_gorilla.txt", "imgt/IGH_human.txt", "imgt/IGH_mouse.txt","imgt/IGH_dog.txt", "imgt/IGH_pig.txt",
             "imgt/TRA_TRD_gorilla.txt", "imgt/TRA_TRD_human.txt", "imgt/TRA_TRD_mouse.txt","imgt/TRA_TRD_dog.txt",
             "imgt/TRB_gorilla.txt", "imgt/TRB_human.txt", "imgt/TRB_mouse.txt","imgt/TRB_dog.txt", "imgt/TRB_pig.txt",]
parse_multiple_embl(file_list, "original")"""

#types = list_feature_types("imgt/IGL_gorilla.txt")
#print(types)
