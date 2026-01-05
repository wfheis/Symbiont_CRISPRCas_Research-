"""
Will Heisler
Script to extract CRISPR array information from a GFF file and write them to a text file.
10/27/2025
"""

import sys
import re

def gff_to_array_info(input_file, output_file, min_spacers=1):
    """
    Extracts full CRISPR array info lines from a combined multi-strain GFF file.
    For each strain, only arrays with more than `min_spacers` spacers are included.

    Output format:
    original header line beginning with ##
    > CRISPR_line
    > CRISPR_line
    """

    strain_arrays = {}        # strain - set of valid array names
    strain_info_lines = {}    # strain - list of CRISPR info lines
    strain_headers = {}       # strain - header line from the GFF file

    # ---------- PASS 1: Identify valid arrays per strain ----------
    strain_name = "unknown_strain"
    with open(input_file, 'r') as f:
        for line in f:
            # Detect strain header
            if line.startswith("##"):
                strain_match = re.search(r"strain[:=]\s*([A-Za-z0-9_-]+)", line)
                if strain_match:
                    strain_name = strain_match.group(1)
                    strain_headers[strain_name] = line.strip()
                    if strain_name not in strain_arrays:
                        strain_arrays[strain_name] = set()
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]

            if feature_type == "CRISPR":
                num_match = re.search(r"Number_of_spacers=(\d+)", attributes)
                name_match = re.search(r"Name=([^;]+)", attributes)
                if num_match and name_match:
                    num_spacers = int(num_match.group(1))
                    array_name = name_match.group(1)
                    if num_spacers > min_spacers: # only consider arrays with more than min_spacers
                        strain_arrays[strain_name].add(array_name)

    # ---------- PASS 2: Collect CRISPR info lines ----------
    current_strain = "unknown_strain"
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("##"): # Detect strain header
                strain_match = re.search(r"strain[:=]\s*([A-Za-z0-9_-]+)", line)
                if strain_match:
                    current_strain = strain_match.group(1)
                    if current_strain not in strain_info_lines: 
                        strain_info_lines[current_strain] = []
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]

            if feature_type == "CRISPR":
                name_match = re.search(r"Name=([^;]+)", attributes)
                if name_match:
                    array_name = name_match.group(1)
                    if current_strain in strain_arrays and array_name in strain_arrays[current_strain]:
                        strain_info_lines[current_strain].append(line.strip())

    # ---------- PASS 3: Write formatted output ----------
    with open(output_file, 'w') as out:
        for strain, crispr_lines in strain_info_lines.items():
            if strain in strain_headers:
                out.write(f"{strain_headers[strain]}\n")
            else:
                out.write(f"##gff-version 3; strain: {strain}\n")

            for crispr in crispr_lines:
                out.write(f"> {crispr}\n")
            out.write("\n")


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python gff_to_array_info.py <input_file> <output_file> [min_spacers]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    min_spacers = int(sys.argv[3]) if len(sys.argv) == 4 else 1

    gff_to_array_info(input_file, output_file, min_spacers)


