"""
Will Heisler
Script to extract spacer sequences from a GFF file and write them to a FASTA file.
10/20/2025
"""
import sys
import re

def gff_to_fasta(input_file, output_file, min_spacers=1):
    """
    Extracts spacer sequences from a combined multi-strain GFF file and writes them to a FASTA file.
    For each strain, only spacers belonging to arrays with more than `min_spacers` spacers are included.
    The FASTA headers include the strain name as a prefix (e.g., >bh11_sp_2039887).
    """

    valid_arrays = set()
    strain_name = "unknown_strain"
    strain_arrays = {}  # maps strain - set of valid arrays

    # ---------- PASS 1: Identify valid arrays per strain ----------
    with open(input_file, 'r') as f:
        for line in f:
            # Detect new strain header line
            strain_match = re.search(r"strain[:=]\s*([A-Za-z0-9_-]+)", line)
            if strain_match:
                strain_name = strain_match.group(1)
                if strain_name not in strain_arrays:
                    strain_arrays[strain_name] = set()
                continue

            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]

            # Identify valid CRISPR arrays
            if feature_type == "CRISPR":
                num_match = re.search(r"Number_of_spacers=(\d+)", attributes) 
                name_match = re.search(r"Name=([^;]+)", attributes)
                if num_match and name_match:
                    num_spacers = int(num_match.group(1))
                    array_name = name_match.group(1)
                    if num_spacers > min_spacers : # SPACER FILTERING CONDITION
                        strain_arrays[strain_name].add(array_name)

    # ---------- PASS 2: Write spacers for valid arrays ----------
    current_strain = "unknown_strain"

    with open(input_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            # Update strain whenever a new header is encountered
            strain_match = re.search(r"strain[:=]\s*([A-Za-z0-9_-]+)", line)
            if strain_match:
                current_strain = strain_match.group(1)
                continue

            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "CRISPRspacer":
                continue

            attributes = fields[8]
            seq_match = re.search(r"sequence=([^;]+)", attributes)
            id_match = re.search(r"ID=([^;]+)", attributes)
            parent_match = re.search(r"Parent=([^;]+)", attributes)

            if seq_match and id_match and parent_match:
                sequence = seq_match.group(1)
                spacer_id = id_match.group(1)
                parent_id = parent_match.group(1)

                # Check if spacer belongs to a valid array for the current strain
                if current_strain in strain_arrays and parent_id in strain_arrays[current_strain]:
                    out.write(f">{current_strain}_{spacer_id}\n{sequence}\n")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python gff_to_spacers_fasta.py <input_file> <output_file> [min_spacers]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    min_spacers = int(sys.argv[3]) if len(sys.argv) == 4 else 1

    gff_to_fasta(input_file, output_file, min_spacers)