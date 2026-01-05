
def reformat_fasta_headers(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Extract the scientific name in brackets
                if '[' in line and ']' in line:
                    name_start = line.index('[') + 1
                    name_end = line.index(']')
                    sci_name = line[name_start:name_end].replace(' ', '_')
                    
                    # Remove original scientific name and reformat the header
                    header_without_brackets = line[:line.index('[')].strip()
                    new_header = f">{sci_name} {header_without_brackets[1:]}\n"
                    outfile.write(new_header)
                else:
                    # Write unmodified if no brackets
                    outfile.write(line)
            else:
                # Write sequence lines as-is
                outfile.write(line)


reformat_fasta_headers('csy4F_prot_large_mafft', 'csy4F_prot_large_mafft_reform')