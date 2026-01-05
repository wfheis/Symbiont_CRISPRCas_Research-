"""
Will Heisler
Script to extract representative DR sequence from a GFF file and write them to a FASTA file.
10/20/2025
"""
import sys
import re

def gff_to_dr_fasta(input_file, output_file, min_spacers=1):
    """
    Extracts DR sequences from 'CRISPR' feature lines in a combined multi-strain GFF file
    and writes them to a FASTA file.

    For each strain, only CRISPR arrays with more than `min_spacers` spacers are included.
    The FASTA headers include the strain name and CRISPR Name field 
    (e.g., >bh11_1_1929629_1930856).
    """

    strain_name = "unknown_strain"

    with open(input_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            # Detect strain header lines (e.g., ##gff-version 3; ... strain: bb410)
            strain_match = re.search(r"strain[:=]\s*([A-Za-z0-9_-]+)", line)
            if strain_match:
                strain_name = strain_match.group(1)
                continue

            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "CRISPR":
                continue

            attributes = fields[8]

            # Extract number of spacers, DR sequence, and CRISPR name
            num_match = re.search(r"Number_of_spacers=(\d+)", attributes)
            dr_match = re.search(r"DR=([^;]+)", attributes)
            name_match = re.search(r"Name=([^;]+)", attributes)

            if not (num_match and dr_match and name_match):
                continue

            num_spacers = int(num_match.group(1))
            dr_seq = dr_match.group(1)
            crispr_name = name_match.group(1)

            # Apply spacer count filter
            if num_spacers > min_spacers:
                header = f">{strain_name}_{crispr_name}"
                out.write(f"{header}\n{dr_seq}\n")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python gff_to_dr_fasta.py <input_file> <output_file> [min_spacers]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    min_spacers = int(sys.argv[3]) if len(sys.argv) == 4 else 1

    gff_to_dr_fasta(input_file, output_file, min_spacers)


