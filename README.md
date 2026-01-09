# Symbiont_CRISPRCas_Research-
This research project was undertaken with the mentorship of Dr. Suegene Noh at Colby College. The goal of this study was to identify evolutionary patterns in the CRISPR-Cas defense systems possessed by Paraburkholderia symbionts of Dictyostelium discoideum.
> These bacteria are facultative symbionts of the amoeba D. discoideum, which is a useful model for animal immune systems. Therefore, studying these symbionts offers an informative window into how pathogen-host interactions shape bacterial genomes and guide evolutionary mechanisms.
> My particular study focuses on the role played by bacterial defense systems in the dynamics and observed effects of this symbiosis. Bacteriophages are a vital driver of bacterial evolution, and the CRISPR-Cas systems meant to defend against phage infection face unique evolutionary pressures in this symbiotic relationship. Therefore, I am seeking to understand what specific pressures are leading to the conservation of CRISPR-Cas in these symbionts and what their presence and makeup can tell us about the evolutionary history of this symbiosis.

### Scripts
#### `fetch_protien_batch.sh` 
> This shell script was written to take a text file input containing NCBI protein accession IDs, fetch the protein sequences using Entrez Direct, and return the sequences in FASTA format.
#### `gff_to_spacers_fasta.py` , `gff_to_DR_fasta.py` , `gff_to_array_info.py`
> These Python scripts were written to extract key data, such as CRISPR array spacer sequences, direct repeat (DR) sequences, and array metadata, from .gff files produced by running the CRISPRCasFinder software on symbiont genomes. Once extracted, this data was used to compare CRISPR spacer content between strains
#### `build_spacer_network.py` 
> This Python script processes the output of an all-against-all nucleotide BLAST of strain-specific CRISPR array spacers and constructs a network visualizing the number of spacers shared between strains using the NetworkX package.
