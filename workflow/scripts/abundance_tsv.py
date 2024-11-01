import sys
from Bio import SeqIO

def abundance_tsv(input_fasta, output_tsv, num_sequences, abundances):
    
    # Check that relative abundance sums to 1
    if sum(abundances) != 1:
        print(f"Error: The relative abundance values must sum to 1. Current sum: {sum(abundances)}.")
        sys.exit(1)

    # Read the sequences from the input fasta file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    if len(sequences) != num_sequences:
        print(f"Error: The input fasta file contains {len(sequences)} sequences, but {num_sequences} sequences were expected.")
        sys.exit(1)

    # Write the SeqID and relative abundance to the output TSV file
    with open(output_tsv, "w") as out_tsv:
        for i, seq_record in enumerate(sequences):
            seq_id = seq_record.id
            abundance = abundances[i]
            out_tsv.write(f"{seq_id}\t{abundance}\n")

    print(f"Successfully wrote abundance TSV to {output_tsv}.")

# Main script logic
if __name__ == "__main__":
    input_fasta = snakemake.input[0]
    output_tsv = snakemake.output[0]
    input_fasta_renamed = snakemake.input[1]
    output_tsv_renamed = snakemake.output[1]
    num_sequences = snakemake.params.hap_count
    abundances = snakemake.params.abundance
    
    # Call the function
    abundance_tsv(input_fasta, output_tsv, num_sequences, abundances)
    abundance_tsv(input_fasta_renamed, output_tsv_renamed, num_sequences, abundances)

